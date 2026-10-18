#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_sarah_cpg_rows(path, stats):
    """
    Reproduce Sarah's CpG-cytosine filtering.

    A row is retained if it is adjacent by exactly 1 bp to the
    previous or next row on the same chromosome.

    The input does NOT need chromosomes in lexical order.
    We require only:
      1. each chromosome/scaffold occurs in one contiguous block;
      2. positions increase within that block.
    """

    raw_rows = 0
    kept_rows = 0

    previous = None
    previous_has_left_pair = False

    current_chrom = None
    seen_chroms = set()
    previous_pos = None

    with open_text(path) as fh:
        for line_no, line in enumerate(fh, start=1):

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 6:
                raise RuntimeError(
                    f"{path}: line {line_no} has fewer than 6 columns"
                )

            chrom = fields[0]

            try:
                pos1 = int(fields[1])
                methC = int(fields[4])
                unmethC = int(fields[5])
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: invalid numeric field at line {line_no}"
                ) from exc

            # Validate chromosome blocks and within-block order.
            if chrom != current_chrom:
                if chrom in seen_chroms:
                    raise RuntimeError(
                        f"{path}: chromosome {chrom} appears in multiple blocks"
                    )

                seen_chroms.add(chrom)
                current_chrom = chrom
                previous_pos = pos1

            else:
                if pos1 < previous_pos:
                    raise RuntimeError(
                        f"{path}: positions decrease within chromosome "
                        f"{chrom} at line {line_no}"
                    )

                previous_pos = pos1

            current = (chrom, pos1, methC, unmethC)
            raw_rows += 1

            if previous is None:
                previous = current
                previous_has_left_pair = False
                continue

            adjacent = (
                current[0] == previous[0]
                and current[1] - previous[1] == 1
            )

            if previous_has_left_pair or adjacent:
                kept_rows += 1
                yield previous

            previous = current
            previous_has_left_pair = adjacent

    if previous is not None and previous_has_left_pair:
        kept_rows += 1
        yield previous

    stats["raw_rows"] = raw_rows
    stats["kept_rows"] = kept_rows


def read_manifest(path):
    populations = defaultdict(list)

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"population", "individual", "cov_path"}

        if set(reader.fieldnames or []) != required:
            raise RuntimeError(
                "Manifest columns must be exactly: "
                "population, individual, cov_path"
            )

        for row in reader:
            cov_path = Path(row["cov_path"])

            if not cov_path.is_file():
                raise FileNotFoundError(cov_path)

            populations[row["population"]].append(
                (row["individual"], cov_path)
            )

    return populations


def pool_population(population, individuals, output_dir, sort_memory="4G"):
    """
    Filter each individual independently, then reproduce R aggregate()
    using GNU sort + streaming reduction.

    This is independent of chromosome-block order in the original files.
    """

    print(f"\n===== {population} =====")

    if shutil.which("sort") is None:
        raise RuntimeError("GNU/coreutils sort was not found in PATH")

    individual_stats = {}

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{population}.pooled_cpg.tsv.gz"

    with tempfile.TemporaryDirectory(
        prefix=f"{population}_sort_",
        dir=output_dir
    ) as tmpdir:

        env = os.environ.copy()
        env["LC_ALL"] = "C"

        sort_cmd = [
            "sort",
            "-T", tmpdir,
            "-S", sort_memory,
            "-k1,1",
            "-k2,2n",
        ]

        print("Starting external sort...")

        proc = subprocess.Popen(
            sort_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
            env=env,
        )

        try:
            for individual, path in individuals:
                print(f"Filtering: {individual} -> {path}")

                stats = {}
                individual_stats[individual] = stats

                for chrom, pos1, methC, unmethC in iter_sarah_cpg_rows(
                    path, stats
                ):
                    proc.stdin.write(
                        f"{chrom}\t{pos1}\t{methC}\t{unmethC}\n"
                    )

            proc.stdin.close()

            pooled_sites = 0

            with gzip.open(output_path, "wt") as out:
                out.write("chrom\tstart\tend\tmethC\tunmethC\n")

                previous_key = None
                sum_meth = 0
                sum_unmeth = 0

                for line in proc.stdout:
                    chrom, pos, methC, unmethC = line.rstrip("\n").split("\t")

                    pos1 = int(pos)
                    key = (chrom, pos1)

                    if previous_key is not None and key != previous_key:
                        prev_chrom, prev_pos1 = previous_key

                        out.write(
                            f"{prev_chrom}\t{prev_pos1 - 1}\t{prev_pos1}"
                            f"\t{sum_meth}\t{sum_unmeth}\n"
                        )

                        pooled_sites += 1

                        if pooled_sites % 5_000_000 == 0:
                            print(
                                f"  pooled sites written: "
                                f"{pooled_sites:,}"
                            )

                        sum_meth = 0
                        sum_unmeth = 0

                    sum_meth += int(methC)
                    sum_unmeth += int(unmethC)
                    previous_key = key

                if previous_key is not None:
                    prev_chrom, prev_pos1 = previous_key

                    out.write(
                        f"{prev_chrom}\t{prev_pos1 - 1}\t{prev_pos1}"
                        f"\t{sum_meth}\t{sum_unmeth}\n"
                    )

                    pooled_sites += 1

            return_code = proc.wait()

            if return_code != 0:
                raise RuntimeError(
                    f"sort failed with exit code {return_code}"
                )

        finally:
            if proc.stdin and not proc.stdin.closed:
                proc.stdin.close()

        print(f"Output: {output_path}")
        print(f"Pooled CpG cytosine positions: {pooled_sites:,}")

        for individual, _ in individuals:
            stats = individual_stats[individual]

            raw = stats["raw_rows"]
            kept = stats["kept_rows"]
            prop = kept / raw if raw else 0

            print(
                f"{individual}: "
                f"raw={raw:,}, "
                f"kept={kept:,}, "
                f"kept_fraction={prop:.4f}"
            )

    return output_path, pooled_sites, individual_stats


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--manifest",
        required=True,
        help="Tab-separated input manifest"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for population-level pooled CpG counts"
    )

    parser.add_argument(
        "--sort-memory",
        default="4G",
        help="Memory limit passed to GNU sort, e.g. 4G"
    )

    args = parser.parse_args()

    populations = read_manifest(args.manifest)

    expected = {"benthic", "littoral", "riverine"}

    if set(populations) != expected:
        raise RuntimeError(
            f"Expected populations {sorted(expected)}, "
            f"found {sorted(populations)}"
        )

    for population in ("benthic", "littoral", "riverine"):
        pool_population(
            population,
            populations[population],
            args.output_dir,
            sort_memory=args.sort_memory
        )


if __name__ == "__main__":
    main()

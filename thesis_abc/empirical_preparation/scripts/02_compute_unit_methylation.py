#!/usr/bin/env python3

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


POPULATIONS = ("benthic", "littoral", "riverine")


def open_text(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def read_rows(path, expected_cols):
    rows = []

    with open_text(path, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")

        for lineno, row in enumerate(reader, start=1):
            if not row or all(x == "" for x in row):
                continue

            if len(row) != expected_cols:
                raise ValueError(
                    f"{path}: expected {expected_cols} columns at line "
                    f"{lineno}, found {len(row)}"
                )

            rows.append(row)

    return rows


def read_cgis(path, sarah_compat=False):
    raw = read_rows(path, 5)

    # Sarah's units_meanMeth.R reads the headerless CGI BED with
    # header=TRUE. Therefore the first real CGI is accidentally consumed
    # as a header. We reproduce that only during validation.
    if sarah_compat:
        if not raw:
            raise ValueError(f"{path}: CGI file is empty")
        raw = raw[1:]

    out = []

    for chrom, start, end, name, cpgcount in raw:
        out.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "CpGcount": int(cpgcount),
            }
        )

    return out


def read_bins(path):
    out = []
    prefix = "CpGcount_"

    for chrom, start, end, label in read_rows(path, 4):
        if not label.startswith(prefix):
            raise ValueError(
                f"{path}: unexpected CpG-bin label {label!r}"
            )

        out.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "CpGcount": int(label[len(prefix):]),
            }
        )

    return out


def read_archipelagos(path):
    out = []

    for chrom, start, end, arch_id, cgi_number in read_rows(path, 5):
        out.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "parent_id": arch_id,
                "parent_n": int(cgi_number),
            }
        )

    return out


def read_interarchipelagos(path):
    out = []

    for chrom, start, end, boundary_arch, bin_number in read_rows(path, 5):
        out.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "parent_id": boundary_arch,
                "parent_n": int(bin_number),
            }
        )

    return out


def attach_parent_units(children, parents, unit_type):
    """
    Reproduce Sarah's bed_intersect + within-structure counter logic.

    Each CGI must overlap exactly one archipelago.
    Each CpG bin must overlap exactly one interarchipelago.
    """

    children_by_chrom = defaultdict(list)
    parents_by_chrom = defaultdict(list)

    for child in children:
        children_by_chrom[child["chrom"]].append(child)

    for parent in parents:
        parents_by_chrom[parent["chrom"]].append(parent)

    for chrom in children_by_chrom:
        children_by_chrom[chrom].sort(
            key=lambda x: (x["start"], x["end"])
        )

    for chrom in parents_by_chrom:
        parents_by_chrom[chrom].sort(
            key=lambda x: (x["start"], x["end"])
        )

    matched = []

    for chrom in sorted(children_by_chrom):
        child_list = children_by_chrom[chrom]
        parent_list = parents_by_chrom.get(chrom, [])

        parent_ptr = 0

        for child in child_list:

            while (
                parent_ptr < len(parent_list)
                and parent_list[parent_ptr]["end"] <= child["start"]
            ):
                parent_ptr += 1

            overlaps = []
            j = parent_ptr

            while (
                j < len(parent_list)
                and parent_list[j]["start"] < child["end"]
            ):
                parent = parent_list[j]

                if parent["end"] > child["start"]:
                    overlaps.append(parent)

                j += 1

            if len(overlaps) != 1:
                raise ValueError(
                    f"{unit_type} unit "
                    f"{child['chrom']}:{child['start']}-{child['end']} "
                    f"overlaps {len(overlaps)} parent structures; "
                    f"expected exactly 1"
                )

            parent = overlaps[0]

            overlap_bp = (
                min(child["end"], parent["end"])
                - max(child["start"], parent["start"])
            )

            matched.append(
                {
                    **child,
                    "parent_id": parent["parent_id"],
                    "parent_start": parent["start"],
                    "parent_end": parent["end"],
                    "parent_n": parent["parent_n"],
                    "bp_number": overlap_bp,
                }
            )

    # Sarah's ave(..., parent_id, FUN=seq_along) counters are generated
    # before archipelago/interarchipelago units are combined.
    matched.sort(
        key=lambda x: (x["chrom"], x["start"], x["end"])
    )

    counters = defaultdict(int)
    out = []

    for item in matched:
        parent_id = item["parent_id"]

        counters[parent_id] += 1
        counter = counters[parent_id]

        if unit_type == "archipelago":
            unit_id = f"{parent_id}_island_{counter}"

            structure_id = (
                f"arch_{parent_id}"
                f"_start_{item['parent_start']}"
                f"_end_{item['parent_end']}"
                f"_islandNumber_{item['parent_n']}"
            )

        elif unit_type == "interarchipelago":
            unit_id = f"{parent_id}_bin_{counter}"

            structure_id = (
                f"interarch_{parent_id}"
                f"_start_{item['parent_start']}"
                f"_end_{item['parent_end']}"
                f"_binNumber_{item['parent_n']}"
            )

        else:
            raise ValueError(f"Unknown unit type: {unit_type}")

        out.append(
            {
                "chrom": item["chrom"],
                "start": item["start"],
                "end": item["end"],
                "bp_number": item["bp_number"],
                "unit_id": unit_id,
                "CpGcount": item["CpGcount"],
                "type": unit_type,
                "structure_id": structure_id,
            }
        )

    return out


def build_unit_skeleton(
    cgi_file,
    bin_file,
    arch_file,
    interarch_file,
    sarah_compat,
):
    cgis = read_cgis(
        cgi_file,
        sarah_compat=sarah_compat,
    )

    bins = read_bins(bin_file)

    arch = read_archipelagos(arch_file)
    interarch = read_interarchipelagos(interarch_file)

    arch_units = attach_parent_units(
        cgis,
        arch,
        "archipelago",
    )

    interarch_units = attach_parent_units(
        bins,
        interarch,
        "interarchipelago",
    )

    # Sarah: rbind(arch_df, interarch_df), then order(chrom, start)
    units = arch_units + interarch_units

    units.sort(
        key=lambda x: (x["chrom"], x["start"], x["end"])
    )

    return units, len(arch_units), len(interarch_units)


def map_population(units, pooled_file):
    """
    Reproduce valr::bed_map using standard BED coordinates.

    The pooled CpG file contains:
      chrom, start, end, methC, unmethC

    where each cytosine is represented as a 1-bp BED interval.
    """

    meth = [0] * len(units)
    unmeth = [0] * len(units)
    hits = [0] * len(units)

    units_by_chrom = defaultdict(list)

    for idx, unit in enumerate(units):
        units_by_chrom[unit["chrom"]].append(idx)

    current_chrom = None
    chrom_units = []
    unit_ptr = 0
    active = []

    previous_chrom = None
    previous_start = None

    processed = 0

    with open_text(pooled_file, "rt") as fh:
        reader = csv.reader(fh, delimiter="\t")

        for lineno, row in enumerate(reader, start=1):

            if not row:
                continue

            # Our Step 1 pooled files contain this header.
            if (
                lineno == 1
                and row[:5]
                == ["chrom", "start", "end", "methC", "unmethC"]
            ):
                continue

            if len(row) != 5:
                raise ValueError(
                    f"{pooled_file}: expected 5 columns at line "
                    f"{lineno}, found {len(row)}"
                )

            chrom, start_s, end_s, meth_s, unmeth_s = row

            start = int(start_s)
            end = int(end_s)
            m = int(meth_s)
            u = int(unmeth_s)

            # Step 1 output is expected to be sorted by chromosome/start.
            if previous_chrom == chrom:
                if start < previous_start:
                    raise ValueError(
                        f"{pooled_file}: positions are not sorted at "
                        f"line {lineno}"
                    )

            elif (
                previous_chrom is not None
                and chrom < previous_chrom
            ):
                raise ValueError(
                    f"{pooled_file}: chromosomes are not lexically "
                    f"sorted at line {lineno}"
                )

            previous_chrom = chrom
            previous_start = start

            if chrom != current_chrom:
                current_chrom = chrom
                chrom_units = units_by_chrom.get(chrom, [])
                unit_ptr = 0
                active = []

            # Remove units ending before this CpG interval.
            active = [
                idx
                for idx in active
                if units[idx]["end"] > start
            ]

            # Add any units that may overlap this CpG interval.
            while unit_ptr < len(chrom_units):
                idx = chrom_units[unit_ptr]

                if units[idx]["start"] >= end:
                    break

                if units[idx]["end"] > start:
                    active.append(idx)

                unit_ptr += 1

            # Count the CpG cytosine in every overlapping unit.
            for idx in active:
                unit = units[idx]

                if (
                    unit["start"] < end
                    and unit["end"] > start
                ):
                    meth[idx] += m
                    unmeth[idx] += u
                    hits[idx] += 1

            processed += 1

            if processed % 5_000_000 == 0:
                print(
                    f"  pooled CpG rows processed: "
                    f"{processed:,}",
                    flush=True,
                )

    return meth, unmeth, hits


def write_population(
    units,
    pooled_file,
    out_file,
):
    meth, unmeth, hits = map_population(
        units,
        pooled_file,
    )

    out_file = Path(out_file)
    out_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    summary = {
        "rows": len(units),
        "non_na": 0,
        "na": 0,
        "arch_non_na": 0,
        "inter_non_na": 0,
    }

    with gzip.open(
        out_file,
        "wt",
        newline="",
    ) as fh:

        writer = csv.writer(
            fh,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "bp_number",
                "unit_id",
                "coverage",
                "CpGcount",
                "type",
                "structure_id",
                "meth_frac",
            ]
        )

        for idx, unit in enumerate(units):

            coverage_value = (
                meth[idx] + unmeth[idx]
            )

            # valr::bed_map returns missing values when the
            # unit has no overlapping pooled CpG observations.
            if hits[idx] == 0 or coverage_value == 0:
                coverage = "NA"
                meth_frac = "NA"

                summary["na"] += 1

            else:
                coverage = coverage_value
                meth_frac = repr(
                    meth[idx] / coverage_value
                )

                summary["non_na"] += 1

                if unit["type"] == "archipelago":
                    summary["arch_non_na"] += 1
                else:
                    summary["inter_non_na"] += 1

            writer.writerow(
                [
                    unit["chrom"],
                    unit["start"],
                    unit["end"],
                    unit["bp_number"],
                    unit["unit_id"],
                    coverage,
                    unit["CpGcount"],
                    unit["type"],
                    unit["structure_id"],
                    meth_frac,
                ]
            )

    return summary


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute CGI/bin-level methylation from pooled CpG "
            "counts, reproducing Sarah's units_meanMeth.R logic."
        )
    )

    parser.add_argument(
        "--benthic",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--littoral",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--riverine",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--cgi-file",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--cpgbin-file",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--archipelago-file",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--interarchipelago-file",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--sarah-compat",
        action="store_true",
        help=(
            "Reproduce Sarah's CGI-reading bug by dropping "
            "the first CGI row. Use for validation only, "
            "not production."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    units, n_arch, n_inter = build_unit_skeleton(
        args.cgi_file,
        args.cpgbin_file,
        args.archipelago_file,
        args.interarchipelago_file,
        args.sarah_compat,
    )

    expected_arch = (
        48722
        if args.sarah_compat
        else 48723
    )

    expected_inter = 156514
    expected_total = expected_arch + expected_inter

    print(
        f"Sarah compatibility mode: "
        f"{args.sarah_compat}"
    )

    print(
        f"Archipelago units: "
        f"{n_arch:,}"
    )

    print(
        f"Interarchipelago units: "
        f"{n_inter:,}"
    )

    print(
        f"Total units: "
        f"{len(units):,}"
    )

    if (
        n_arch != expected_arch
        or n_inter != expected_inter
        or len(units) != expected_total
    ):
        raise ValueError(
            "Unexpected unit counts: "
            f"arch={n_arch}, "
            f"interarch={n_inter}, "
            f"total={len(units)}; "
            f"expected arch={expected_arch}, "
            f"interarch={expected_inter}, "
            f"total={expected_total}"
        )

    pooled = {
        "benthic": args.benthic,
        "littoral": args.littoral,
        "riverine": args.riverine,
    }

    for population in POPULATIONS:

        print(
            f"\n===== {population} =====",
            flush=True,
        )

        out_file = (
            args.output_dir
            / f"{population}.mean_meth_units.tsv.gz"
        )

        summary = write_population(
            units,
            pooled[population],
            out_file,
        )

        print(
            f"Output: {out_file}"
        )

        print(
            f"rows: "
            f"{summary['rows']:,}"
        )

        print(
            f"non-NA meth_frac: "
            f"{summary['non_na']:,}"
        )

        print(
            f"NA meth_frac: "
            f"{summary['na']:,}"
        )

        print(
            f"archipelago non-NA: "
            f"{summary['arch_non_na']:,}"
        )

        print(
            f"interarchipelago non-NA: "
            f"{summary['inter_non_na']:,}"
        )


if __name__ == "__main__":
    main()

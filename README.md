This repository contains a Snakemake pipeline for preprocessing whole-genome bisulfite sequencing (WGBS) data.

The pipeline includes:
- Data download from SRA
- Read trimming (TrimGalore)
- Alignment (Bismark)
- Deduplication
- Methylation extraction
- Generation of per-site and 50 bp window methylation levels

The workflow reproduces the preprocessing steps described in Vernaz et al. (2022).

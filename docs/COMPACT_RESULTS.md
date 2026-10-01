# Lossless Compact Results

The upload copy stores CSV files of at least 512 KiB as gzip-compressed
`.csv.gz` files. The initial compaction replaced 13 large CSV files; their
original sizes, compressed sizes and uncompressed MD5 hashes are recorded in
`RESULT_COMPRESSION.csv` in this directory. No observations, attempted fits,
failure records, replications or numeric digits were removed or rounded.

The `results` directory is approximately 17.1 MiB after compaction, instead of
59.4 MiB. The largest individual compressed file is approximately 3.4 MiB.
These sizes describe the prepared snapshot, not a guarantee about a particular
browser upload limit or the size of an existing Git history.

## Reading Results

Every CSV-reading entry point sources `R/results_io.R`. `read_jkss_csv()` reads
the named `.csv` if present, otherwise it reads the `.csv.gz` counterpart through
an R gzip connection. The logical `.csv` names used in the manuscript and other
documentation remain unchanged. Existing `make` and `Rscript` commands work
without manual extraction or additional R packages.

For an individual file:

```r
source("R/results_io.R")
records <- read_jkss_csv("results/main_simulation_records.csv")
```

## After Rerunning Analyses

Analysis scripts continue writing ordinary CSV files. A new plain CSV takes
precedence over an older compressed file. Run `make compact` (or
`Rscript R/compress_results.R`) to replace large CSVs with updated compressed
versions before uploading. The complete `run_all.R` workflow also performs
this compaction before building the reproducibility archive.

Compression uses binary streams and checks restored-file MD5 equality before
removing a redundant plain CSV. The source is rechecked before removal. Small
CSV files, figures, tables, PDFs and RDS files are left unchanged. Compaction is
idempotent when no large plain CSVs remain.

## Updating an Existing Repository

Upload the updated `R/`, `Makefile` and `run_all.R` together with `results/`.
If an uncompressed counterpart from `RESULT_COMPRESSION.csv` is already in the
remote repository, replace it rather than retaining both copies. Do not remove
other summary CSVs. The preparation task does not rewrite Git history, change
remote files or remove the original uncompressed `Latex/JKSS` study records.

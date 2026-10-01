# JKSS Upload Snapshot

Prepared: 2026-09-29; document snapshot refreshed: 2026-09-30.
Source: the author's `Latex/JKSS` project, with manuscript and supplement last
revised on 2026-09-30.

## Included

- The current manuscript, supplement and all required LaTeX modules.
- Original R estimation, simulation, application, plotting and validation code.
- The JOBS II illustration data, complete retained results, seeds, failure
  records, analysis session records and vector figures.
- The checked direct-LaTeX PDFs in `output/pdf/`.
- Study-design records for the application, variance and confounding
  sensitivity experiments, curve simulation and calibrated power comparison.

## Packaging-Only Changes

The statistical routines, study data, manuscript text and stored numerical
values are preserved without scientific changes. Large result CSVs are now
losslessly compressed; CSV-reading entry points use a shared plain/gzip reader.
Existing validation reports may be refreshed by running the validation commands.

- `README.md` describes this JKSS snapshot and the current section numbering.
- `R/check_dependencies.R` checks the required R packages without installing them.
- `Makefile` adds a dependency check and includes effect-map and application
  checks in `make validate`.
- `run_all.R` uses a portable repository-root message, checks R dependencies,
  includes the new documentation/PDFs in the reproducibility ZIP, and copies
  direct LaTeX PDF output instead of rewriting it through Ghostscript. The
  latter avoids the font corruption observed in the previous postprocessing.
- Citation metadata, upload instructions, Git exclusions and data notices were
  added. No publication status or DOI is asserted.
- The compact revision adds `R/results_io.R`, `R/compress_results.R` and
  `R/test_results_io.R`. Eighteen CSV-reading scripts use the shared reader;
  result-presence checks also recognize compressed counterparts. Estimation
  equations, random-number generation and analysis formulas are unchanged.
  `make compact` and the final archive step perform verified lossless compression.

## Feasible-Root Theory Update

The 2026-09-29 revision connects the implemented Gaussian-kernel estimator to
Kim (2023, Theorem 2), proves first-order equivalence of the full-kernel and
leave-one-out equations under stated central-set and tail conditions, and gives
sufficient conditions for the deterministic rule to select the consistent
root. The density floor, solver tolerance, score-acceptance threshold and
finite-difference step are now implemented as capped asymptotic sequences. At
every sample size used in the article, their values equal the previous constants,
so the retained estimates, tables and figures are unchanged.

## Supplement Order and Interpretation Update

The 2026-09-30 document revision aligns the Supplementary Material with the main article:
theory, main/comparator simulation and paired audit, calibrated power and
numerical sensitivity, variance misspecification, confounding-sensitivity
methodology and simulation, JOBS II, then reproducibility. All 16 tables and six
figures are retained. The main text and generated table notes use the revised
numbers. Added explanations distinguish model-based precision and weak-effect
power gains from numerical reliability, pointwise calibration and empirical
diagnostics. Stored estimates, data-generating mechanisms and Monte Carlo
replication counts are unchanged.

## Excluded

Cover letters, referee reports, internal journal-selection notes, earlier
drafts, alternative heteroscedastic projects, build logs, local R libraries,
temporary checkpoints and duplicate top-level PDFs are not part of the upload
snapshot. The separately maintained R package is also not bundled: all code
needed for this manuscript is self-contained in `R/`.

The previous `github` directory remains unchanged. The 2026-09-30 document
revision is synchronized with the source JKSS project. No remote repository
update or release was performed by this local revision.

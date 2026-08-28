# Semiparametrically Efficient Causal Mediation

This repository contains the complete reproducibility materials for
**Semiparametrically Efficient Inference for Causal Mediation Effects in Linear
Models with Unspecified Error Distributions**.

**Author:** Mijeong Kim  
**Affiliation:** Department of Statistics, Ewha Womans University, Seoul,
Republic of Korea  
**ORCID:** [0000-0002-3578-5413](https://orcid.org/0000-0002-3578-5413)  
**Contact:** [m.kim@ewha.ac.kr](mailto:m.kim@ewha.ac.kr)

The repository includes executable R code, fixed random-number seeds,
replication-level Monte Carlo records, numerical-failure logs, generated tables
and figures, manuscript sources, package versions, and session information.
Numerical failures are retained rather than silently replaced by OLS fits.

## Repository structure

- `R/semiparametric_mediation.R`: efficient scores, estimators, natural-effect
  map, data generators, and comparator wrappers
- `R/run_simulations.R`: main, comparator, size-power, and algorithm-sensitivity
  experiments
- `R/run_application.R`: JOBS II analysis and residual diagnostics
- `R/validate_outputs.R`: record-count, metric-range, and effect-decomposition
  checks
- `R/make_manuscript_assets.R`: deterministic table and figure generation
- `R/smoke_test.R`: small end-to-end estimator check
- `results/`: full retained numerical records, summaries, status logs, generated
  LaTeX tables, package versions, and session information
- `figures/`: vector figures used in the article and supplement
- `manuscript.tex` and `supplement.tex`: article and supplementary sources
- `run_all.R`: complete one-command workflow

## Requirements

- R 4.3 or later
- R packages: `Matrix`, `rootSolve`, `sn`, `mediation`, and `medflex`
- A LaTeX installation with `latexmk` or `pdflatex` to rebuild the PDFs

Install the non-recommended R packages with:

```r
install.packages(c("rootSolve", "sn", "mediation", "medflex"))
```

## Quick verification

Run the small estimator-level check:

```sh
Rscript R/smoke_test.R
```

Validate the supplied full numerical output:

```sh
Rscript R/validate_outputs.R
```

The validation checks 16,000 main-simulation status records, 10,000 comparator
status records, reported metric ranges, application output, and both total-effect
decompositions.

## Rebuild from supplied results

To regenerate the application output, all tables and figures, the article and
supplement PDFs, and the reproducibility archive without rerunning the Monte
Carlo experiments, run:

```sh
CS_SKIP_SIMULATIONS=1 Rscript run_all.R
```

## Full reproduction

Run the complete workflow from this repository root:

```sh
Rscript run_all.R
```

This command:

1. Runs all Monte Carlo experiments with the reported replication counts.
2. Reanalyzes the JOBS II data distributed with `mediation`.
3. Validates all retained output.
4. Regenerates every table and figure.
5. Compiles the article and supplement when LaTeX is available.
6. Creates `output/CS_Online_Resource_2_reproducibility.zip`.

The comparator experiment is computationally intensive: it analyzes 2,000 data
sets and uses 499 nonparametric bootstrap resamples per data set. Independent
Monte Carlo tasks are parallelized on Unix-like systems.

For a short development run, override the replication counts and worker count:

```sh
CS_REPS_MAIN=2 \
CS_REPS_COMPARATOR=2 \
CS_REPS_POWER=2 \
CS_REPS_SENSITIVITY=2 \
CS_CORES=2 \
Rscript R/run_simulations.R
```

This command overwrites the supplied full numerical CSV files. Use it only in a
working copy or restore the full output before rebuilding manuscript assets.

## Random-number seeds

- Main simulation: `20260827`
- Comparator simulation: `20260828`
- Size and power simulation: `20260829`
- Algorithm sensitivity: `20260830`
- JOBS II quasi-Bayesian analysis: `20260831`
- JOBS II bootstrap analysis: `20260901`
- Error-distribution figure: `20260902`

Each parallel simulation task receives a stored task-specific seed, so simulated
data do not depend on the number of forked workers.

## Data

The JOBS II data are distributed with the R package `mediation`. The application
script loads the data directly from that package; no restricted or identifiable
data are included in this repository.

## Citation

Please cite the accompanying manuscript. Publication metadata will be added
after acceptance.

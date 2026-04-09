# Semiparametric Causal Mediation Analysis for Linear Models with Non-Gaussian Errors

This repository contains the code, derived results, and manuscript files for the Statistics in Medicine submission

`Semiparametric Causal Mediation Analysis for Linear Models with Non-Gaussian Errors: Applications to Drug Treatment and Social Program Evaluation`.

The project studies semiparametric causal mediation analysis for linear models with possibly non-Gaussian errors, with applications to the `uis` drug-treatment data and the `jobs` social-program data.

## Repository layout

- `scripts/stat_med_mediation_revision.R`: main analysis script for simulations, the near-boundary power example, and the two empirical applications
- `scripts/stat_med_figures.R`: figure-generation script for the application-effect figure and supplementary residual histograms
- `results/`: derived CSV outputs used in the manuscript
- `figures/`: figure files used in the manuscript and supplement
- `stat_med_submission.tex`: main manuscript
- `stat_med_submission_supplement.tex`: online supplementary material
- `stat_med_coverletter.tex`: cover letter

## Reproducibility

The analysis scripts rely on the supporting code in `R/jobs_semipara_delta_revision.R` together with several R packages, including `Matrix`, `sn`, `quantreg`, and `mediation`.

Running `scripts/stat_med_mediation_revision.R` regenerates the main result tables saved under `results/`. Running `scripts/stat_med_figures.R` regenerates the figure PDFs saved under `figures/`.

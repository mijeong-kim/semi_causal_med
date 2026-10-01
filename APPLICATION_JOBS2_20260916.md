# JOBS II application record

## Data provenance

The application uses the `jobs` data object distributed with version 4.5.1 of
the R package `mediation`. The exact object used for the manuscript is stored in
`data/jobs.RData`. It contains 899 observations and no missing values among the
variables selected for analysis.

## Analysis definition

- Treatment: randomized intervention assignment (`treat`).
- Mediator: job-search self-efficacy (`job_seek`).
- Outcome: post-intervention depressive symptoms (`depress2`).
- Baseline adjustment: pretreatment depression (`depress1`), economic hardship
  (`econ_hard`), sex and age.
- Outcome mean: treatment, mediator, treatment-by-mediator interaction and all
  baseline terms.
- Target population: the empirical baseline-covariate distribution of the 899
  package observations.

The application compares OLS delta inference, the density-score proposal, the
quasi-Bayesian and percentile-bootstrap options in `mediation`, and the natural
effect model in `medflex`. Seeds are 20260831 for quasi-Bayesian simulation and
20260901 for the bootstrap.

## Prespecified checks

Residual density and variance diagnostics are reported for both regressions.
The covariate-sensitivity analysis removes pretreatment depression while
retaining the sample, other regressors and numerical settings. The algorithm
check uses starting-value scales 0.05, 0.10 and 0.20.

The reduced-form confounding analysis uses the same observations and baseline
variables, separate outcome regressions by treatment stratum, a common
analyst-specified correlation over [-0.8, 0.8], and both signs of the residual
R-squared calibration.

## Interpretation

The data provide a reproducible behavioral-intervention illustration. Residual
diagnostics show both non-Gaussian shape and variance heterogeneity, so nominal
intervals are interpreted under the stated working model and alongside the
variance-misspecification experiment. Causal interpretation additionally relies
on the identification and sensitivity assumptions stated in the manuscript.

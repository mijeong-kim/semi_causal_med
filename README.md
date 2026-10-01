# Semiparametric Causal Mediation: Code and Reproducibility Materials

**Target journal:** Journal of the Korean Statistical Society  
**Author:** Mijeong Kim  
**Affiliation:** Department of Statistics, Ewha Womans University, Seoul, Republic of Korea  
**Correspondence:** m.kim@ewha.ac.kr  
**ORCID:** 0000-0002-3578-5413  
**Public repository:** https://github.com/mijeong-kim/semi_causal_med

This repository provides analysis code, data, retained numerical results and
standalone figures for **Semiparametric Inference and Sensitivity Analysis
for Causal Mediation in Linear Models with Unspecified Error Distributions**.
The replication snapshot was refreshed on **2026-10-01**. It concerns the
independent-error, density-score mediation study, not the separate cross-fitted
heteroscedastic project.

## Public Repository Scope

The main article and its Supplementary Material are reserved for separate
submission to the Journal of the Korean Statistical Society. Their PDFs and
LaTeX document sources are not part of the public distribution. This repository
is for computational reproducibility, not for distributing the submission
documents.

The public materials are the R code in `R/`, study data in `data/`, numerical
records and diagnostics in `results/`, standalone figures in `figures/`, and
the instructions and provenance needed to use them. Figure PDFs are individual
plots, not the article or supplementary document. The author-side `output/`
directory, including compiled documents and backup ZIPs, is not an upload target.

The commands below reproduce analyses, numerical summaries and figures without
the article or supplementary sources. They also generate table fragments locally
in `results/`; these generated `.tex` fragments are numerical exports, not the
LaTeX sources of the submission documents.

**Compact results:** large CSV files are supplied as lossless `.csv.gz` files.
All replication rows and original numeric precision are retained. The supplied
R scripts read these files directly without extracting them.
Numerical summaries remain available as plain CSV files
except for the large full-curve summary. See [compression details](docs/COMPACT_RESULTS.md).

- [Standalone figures](figures/)
- [Numerical results and diagnostics](results/)
- [Executed checks and their scope](docs/VERIFICATION.md)
- [Data provenance](data/README.md) and [third-party notices](THIRD_PARTY_NOTICES.md)

The main experiments compare OLS and the proposed estimator at n=200 and 500
with 1,000 replications per error law. The separate n=300 comparator experiment
uses 1,000 replications per law, 1,000 quasi-Bayesian draws and 499 bootstrap
resamples. The materials also include independently calibrated power comparisons
with Huber-based estimators, variance-misspecification experiments, confounding
sensitivity curves, and the JOBS II illustration. Failures and seed records are
retained alongside successful-fit summaries.

The implementation combines regression density scores, stacked inference and
a deterministic multi-start algorithm. The sensitivity extension retains
Imai et al.'s fixed-correlation identification map while changing nuisance
estimation and propagating its joint covariance. The scope and assumptions of
each experiment are summarized below; the full methodological exposition is
part of the separately submitted documents.

## Quick Start

Run the following from the repository root, after installing the requirements:

```sh
Rscript R/check_dependencies.R
make smoke
make validate
make assets
```

`make smoke` runs small estimator and code checks, and `make validate` checks
the retained numerical results. `make assets` rebuilds tables and figures from
those results without rerunning the full Monte Carlo experiments or the
application. These commands do not compile an article or create a ZIP and do
not need the top-level `.tex` files. Without `make`, run the corresponding
`Rscript` commands in `Makefile`.

Use the explicit targets above rather than bare `make`, `make all`, `make pdf`
or `Rscript run_all.R`. Those legacy author-side entry points compile the
article and supplement and require the complete LaTeX source tree; they are
not the public analysis-only workflow.

## Requirements

- R 4.3 or later
- R packages: `Matrix`, `rootSolve`, `sn`, `mediation`, `medflex`, and `quantreg`
- `make` for the convenience targets, or run their `Rscript` commands directly

LaTeX and `zip` are not required for the analysis-only workflow. Figure PDFs
are generated directly by R; table `.tex` files can be generated without
compiling an article.

Install the non-recommended R packages with:

```r
install.packages(c("rootSolve", "sn", "mediation", "medflex", "quantreg"))
```

`Matrix` is normally included with R; install a version compatible with your R
installation if the dependency check reports it missing. The scripts check
dependencies but never install packages automatically. Recorded analysis
versions appear in `results/package_versions.csv` and the experiment-specific
`*_sessionInfo.txt` files; they are provenance records, not a dependency lockfile.

## Reproduce the Analyses

To rerun all Monte Carlo experiments and both JOBS II analyses, use the following
commands in order from the repository root. This replaces the retained numerical
outputs and can take substantial time; use Quick Start to reuse them instead.

```sh
make simulations calibrated-power variance confounding curves
Rscript R/run_application.R
make validate
make assets
make compact
```

The named targets run the main, comparator, nominal-power, numerical-sensitivity,
calibrated-power, variance-misspecification and confounding-sensitivity studies.
The `confounding` target also runs the reduced-form JOBS II sensitivity analysis;
`R/run_application.R` runs the primary JOBS II comparison. Validation checks the
effect-map derivatives, replication counts, summaries and decompositions.
The final targets regenerate the tables and figures and losslessly compact the
large result files. They do not compile either submission document.
Do not add `make -j`: the experiments have their own worker controls, and the
listed targets must run sequentially because later steps use earlier results.

The comparator experiment is computationally intensive because each of 4,000
data sets uses 499 nonparametric bootstrap resamples. Parallel workers are used
on Unix-like systems.

To run only the comparator experiment at the manuscript settings:

```sh
JKSS_CORES=8 Rscript R/run_comparator.R
```

This runs 1,000 replications per distribution at n=300, with 1,000 quasi-Bayesian
draws and 499 bootstrap resamples within each data set. Source- and seed-checked
checkpoints in `tmp/comparator_checkpoints_1000` permit resumption. Task seeds,
analysis provenance and independent validation are recorded in
`results/comparator_seeds.csv`, `results/comparator_sessionInfo.txt` and
`results/comparator_validation.txt`. Then use `make validate` and `make assets`
to check the results and regenerate tables and figures.

When extending a retained shorter run, `JKSS_COMPARATOR_REUSE_PREFIX=1` checks
every original OLS fit and all five methods at the first, middle and last
replications of each distribution before reusing that prefix. The default
standalone command computes all replications, reusing only matching checkpoints.

To regenerate tables and figures from the supplied full simulation output without
rerunning either the Monte Carlo experiments or the application:

```sh
make assets
```

This target preserves the original numerical records, package-version record and
analysis session information. It regenerates the four-panel error illustration
with its stated seed, without rerunning any mediation experiment. Its environment is
recorded separately in `results/asset_sessionInfo.txt`. Application refits
write their own `results/application_sessionInfo.txt`.

## Independently Calibrated Power and Huber Comparators

The additional 2026-09-23 experiment preserves the original size-power files.
It uses the same asymmetric-mixture design (n=300, gamma=-0.26, eta=0.8),
5,000 independent calibration-null data sets, and 1,000 new evaluation data
sets at each beta2=0,0.1,0.2,0.3,0.4. OLS, Huber-FIX, Huber-SEL and the
unchanged proposed estimator analyze the same task data. Huber comparators
target PNIE using slopes only; their intercepts need not equal mean intercepts.
The tuning rules follow Wang, Peng and Tong (2025), with final-residual
centered-design influence covariance adapted to interaction and joint PNIE
inference. This is not an exact replication of their simple-mediation Sobel code.

The design and tuning rules are implemented in `R/calibrated_power.R` and
`R/run_calibrated_power.R`. Method-specific 95th-percentile cutoffs are computed
only from valid calibration-null fits; held-out null data evaluate the achieved
size. These are DGP-specific simulation
diagnostics, not proposed critical values for applications. Independent resampling
of both Monte Carlo stages preserves method pairing and quantifies calibration
and evaluation uncertainty. Conditional, all-attempt and common-valid results
are retained, including every failure and its message.

```sh
JKSS_CORES=8 make calibrated-power
make assets
```

Master seed: `20260923`; two-stage bootstrap seed: `20260924` (1,000 resamples).
Source-, seed- and package-version-checked checkpoints are stored in
`tmp/power_checkpoints_5000_1000`. The complete records and validation audit are
`results/calibrated_power_records.csv` and `results/calibrated_power_validation.txt`.
The generator now checks for `sn` only when requesting skew-normal errors; no
data-generating formula or random draw order changed. The new power experiment
therefore does not require `sn`, but the complete existing asset workflow does.

Development runs cannot overwrite the manuscript's new power results:

```sh
JKSS_POWER_CALIBRATION=2 JKSS_POWER_EVALUATION=1 \
JKSS_POWER_OUTPUT=tmp/power_smoke JKSS_CORES=2 \
Rscript R/run_calibrated_power.R
```

## Application scope

The application uses the 899-observation JOBS II illustration data distributed
with version 4.5.1 of the `mediation` R package. Treatment is randomized
intervention assignment, the mediator is job-search self-efficacy, and the
outcome is post-intervention depressive symptoms. Both models adjust for
pretreatment depression, economic hardship, sex and age; the outcome model
retains treatment-mediator interaction. Every selected variable is complete.

The data object is included as `data/jobs.RData`, so application reproduction is
offline and version-stable. Residual diagnostics show strong non-Gaussian shape
and evidence of variance heterogeneity. The main analysis therefore reports
working-model intervals together with these diagnostics. A covariate-sensitivity
analysis removes pretreatment depression while retaining the same observations
and remaining model terms. See `APPLICATION_JOBS2_20260916.md` for the analysis
record and provenance.

The imputation-based `medflex` comparator uses `Y ~ T0 * T1 + T0 * X` in the
simulation and the analogous interaction with every centered baseline covariate
in the application. Outcome imputation also includes `T:X` (and all analogous
baseline terms) to satisfy the package's working-model compatibility check;
the extra coefficients are zero in the simulation truth. Its package-native covariance treats the covariate reference
as fixed. The two `mediation` entries are model-based quasi-Bayesian/bootstrap
options, not the separate model-based/design-based frameworks.

## Verification

Run a small estimator-level check with:

```sh
Rscript R/smoke_test.R
Rscript R/test_effect_map.R
Rscript R/test_application.R
Rscript R/validate_outputs.R
Rscript R/test_comparator_outputs.R
Rscript R/paired_fit_audit.R
Rscript R/test_confounding_sensitivity.R
Rscript R/test_confounding_outputs.R
```

The smoke test verifies that OLS, the proposed estimator, both the
quasi-Bayesian and bootstrap `mediation` implementations, and the `medflex`
natural effect model all return the same ordered set of five finite effect
estimates and intervals. The derivative test independently checks propagation
of the estimated baseline-covariate mean and the two exact linear constraints
on the five-effect covariance matrix. The paired-fit audit uses existing
replication identifiers and never creates new data.

The full replication counts and worker count can be overridden for development:

```sh
JKSS_REPS_MAIN=2 \
JKSS_REPS_COMPARATOR=2 \
JKSS_REPS_POWER=2 \
JKSS_REPS_SENSITIVITY=2 \
JKSS_CORES=2 \
Rscript R/run_simulations.R
```

This development command overwrites the full numerical CSV files. Rerun the
default command before validating or regenerating final tables and figures if
reduced counts are used.

## Directory map

- `Makefile`: explicit analysis, validation and table-and-figure targets
- `R/semiparametric_mediation.R`: efficient scores, estimators, effect map, data generators, and comparator wrappers
- `R/run_simulations.R`: main, comparator, size-power, and algorithm-sensitivity experiments
- `R/run_comparator.R`: standalone comparator experiment, resumable checkpoints and validated-prefix extension
- `R/test_comparator_outputs.R`: independent checks of comparator seeds, summaries and reconstructed OLS fits
- `data/jobs.RData`: JOBS II illustration data from `mediation` 4.5.1
- `R/run_application.R`: JOBS II analysis, diagnostics and covariate sensitivity
- `R/make_manuscript_assets.R`: generated LaTeX tables and figures
- `R/paired_fit_audit.R`: matched-data-set comparisons and report-and-cover fractions
- `R/test_effect_map.R`: independent derivative and covariance-constraint tests
- `R/test_application.R`: independent application-model, variance-diagnostic and support checks
- `R/refresh_medflex.R`: version audit refitting only the corrected comparator on original seeds
- `results/`: retained replication-level records, status logs, all reported summaries, generated tables, and `sessionInfo.txt`
- `R/results_io.R`: transparent reading of plain CSV and compressed CSV results
- `R/compress_results.R`: verified lossless compaction after rerunning analyses (`make compact`)
- `figures/`: standalone vector PDF plots; the public graphical outputs
- `docs/`: computational checks, result-compression details and provenance records

The public directory map excludes article and supplementary-document PDFs,
their top-level `.tex` sources, and the author-side `output/` directory. Numerical
exports generated in `results/` are distinct from those submission documents.

## Random-number seeds

- Main simulation: `20260827`
- Comparator simulation: `20260828`
- Size and power simulation: `20260829`
- Algorithm sensitivity: `20260830`
- JOBS II quasi-Bayesian analysis: `20260831`
- JOBS II bootstrap analysis: `20260901`
- Error-distribution figure: `20260902`
- Variance-misspecification sensitivity: `20260917`
- Reduced-form confounding-sensitivity simulation: `20260918`

Each parallel simulation task receives a stored task-specific seed before
execution, making the simulated data invariant to the number of forked workers.

## Variance-misspecification experiment

The design was written in `VARIANCE_SENSITIVITY_PLAN_20260916.md` before examining
the new results. It uses n=300 and 500 attempts per configuration: Gaussian and
asymmetric-mixture innovations, a common constant-variance reference, and
covariate- or treatment-dependent variances at kappa=0.15, 0.30, 0.60. Population
normalization keeps each marginal error variance at one. Mean models, true
natural effects and sequential ignorability are unchanged. The independent-error
restriction is violated intentionally; this is not a new heteroscedastic estimator.

Run only this addition, without refitting earlier studies or the application:

```sh
Rscript R/run_variance_sensitivity.R
Rscript R/test_variance_sensitivity.R
```

`JKSS_VARIANCE_CORES` controls workers (default 4). Checksummed checkpoints in
`tmp/variance_checkpoints_500` permit resumption; delete that checkpoint directory
to force a fresh run. Each of the 1,000 base-data tasks supplies seven scale
configurations with shared random innovations. The output contains 14,000 method
attempts, including failures, and 70,000 five-effect rows. Comparators are the
unchanged proposed method and OLS with the existing stacked empirical sandwich.
Full summaries, paired-fit summaries, seeds, failure messages, MCSEs, a separate
session record, source signatures and verification results are retained in
`results/variance_sensitivity_*`. Earlier results are not reused as new reference
replications. All five effects are available, not just the PNIE/TE table.

For a small isolated development run, use both overrides so the full results
are not overwritten:

```sh
JKSS_VARIANCE_REPS=3 JKSS_VARIANCE_OUT=tmp/variance_smoke Rscript R/run_variance_sensitivity.R
JKSS_VARIANCE_REPS=3 JKSS_VARIANCE_OUT=tmp/variance_smoke Rscript R/test_variance_sensitivity.R
```

The full-run commands above include this experiment; `make validate` checks
its retained outputs and `make assets` rebuilds its tables and figures.
Successful-fit coverage and precision, all-attempted success/report-and-cover,
and paired comparisons must
be interpreted jointly. The experiment supplies the calibration used to interpret
the JOBS II residual-variance diagnostics.
The retained source signature identifies the Monte Carlo run. A subsequent
layout-only adjustment to the coverage figure's outer margin does not change
the generator, estimator or retained numerical results; the complete current
source is supplied for reproduction.

## Reduced-form confounding sensitivity

The design and its pre-simulation pooled-mediator amendment are documented in
`CONFOUNDING_SENSITIVITY_PLAN_20260916.md`. The extension fits pooled
`M ~ T + baseline` and separate `Y ~ baseline` regressions within treatment
strata. It never treats a structural outcome error correlated with the mediator
as independent of that mediator. OLS-RF and DS-RF use the same reduced-form
working model. Only the marginal regression estimating equations differ.
Residual variances, cross-covariances and the empirical baseline distribution
enter the joint sandwich. Imai, Keele and Yamamoto's (2010) fixed-rho map is
unchanged; the code supports separate stratum-specific rho values.

```sh
Rscript R/run_confounding_simulation.R
Rscript R/summarize_confounding.R
Rscript R/run_confounding_application.R
Rscript R/test_confounding_sensitivity.R
Rscript R/test_confounding_outputs.R
Rscript R/make_confounding_assets.R
```

The full simulation has n=300 and 500 attempts for each of two innovation laws
and three true correlations (-0.3, 0, 0.3): 3,000 data sets and 6,000 method fits.
Each fit is evaluated at the true rho and at zero, giving 60,000 effect rows and
12,000 zero-boundary rows. All failures, seeds, summaries, Monte Carlo standard
errors, paired comparisons and session information are in `results/confounding_*`.
Source-signed checkpoints are stored in `tmp/confounding_checkpoints_500`.
`JKSS_CONFOUNDING_CORES` controls workers. For development, set both
`JKSS_CONFOUNDING_REPS=2` and `JKSS_CONFOUNDING_OUT=tmp/confounding_smoke` on the
simulation, summary and output-test commands to preserve the full records.
The application and asset commands always target the manuscript's `results/`.

The JOBS II extension retains exactly the same 899 observations and four baseline
covariates as the primary application. It reports 81 common-rho
values and both signs of a 41-by-41 residual-R-squared grid. The contour axes are
fractions of residual variance (star), not fractions of total variance (tilde).
Their interpretation requires the additional shared-confounder moment
decomposition. Intervals and interval-zero contours are pointwise. The zero of
an estimated effect is not the point at which an interval starts including zero.

Do not interpret the extension as repairing arbitrary heteroscedasticity,
post-treatment confounding or missing outcomes. Common mediator errors and
within-stratum reduced-form error restrictions remain. The zero-boundary
estimators have the same first-order residual-moment representation; improved
curve precision does not imply a more robust causal conclusion. The mixture
simulation shows correlation-dependent gains, not uniform gains; the two
empirical zero boundaries are almost unchanged. The separate variance experiment
also shows severe undercoverage under treatment-dependent scaling, including
some mild departures. Those unfavorable results are retained.

## Expected outputs

### Full-curve simulation extension

The title revision accompanies a full-curve evaluation documented in
`CURVE_SIMULATION_PLAN_20260916.md`. This reuses the SAME 3,000 confounding data
sets (500 per configuration), rather than adding independent replications.
All five effects are evaluated at 33 assumed correlations from -0.8 to 0.8.
Curve bias/RMSE/coverage target the population sensitivity functional at each
assumed correlation. `CausalBias` and `CausalInclusion` instead compare with the
generating causal effect; the targets coincide only at the generating correlation.

```sh
Rscript R/run_confounding_curve.R
Rscript R/summarize_confounding_curve.R
Rscript R/test_confounding_curve.R
Rscript R/test_confounding_curve_outputs.R
Rscript R/make_confounding_curve_assets.R
```

`JKSS_CURVE_CORES` controls workers. Source-signed checkpoints in
`tmp/curve_checkpoints_500` support resumption. For an isolated smoke test, set
both `JKSS_CURVE_REPS=2` and `JKSS_CURVE_OUT=tmp/curve_smoke` for the runner,
summary and output test. For smoke figures also set
`JKSS_CURVE_FIGURES=tmp/curve_smoke/figures` on the asset command.
The output validator compares the refits with the retained original
`results/confounding_*` records and requires those files.

`confounding_curve_basis.csv` stores 30,000 effect-level rows with A, B, VA, CAB
and VB. These exactly reconstruct each estimate as A+h*B and its variance as
VA+2*h*CAB+h^2*VB, where h=rho/sqrt(1-rho^2). Failed attempts remain in the file.
This avoids storing 990,000 redundant grid rows while retaining full numerical
reproducibility. The summary has 1,980 method/effect/grid cells; paired output
has 990 cells. Zero-boundary summaries are also retained. Refits are checked
against every earlier estimate, standard error, status and zero boundary.

The main simulation figure shows population curves, Monte Carlo mean curves
and central 90% empirical sampling ranges. Its shading is NOT a confidence band.
Supplementary figures show pointwise interval coverage and paired RMSE ratios;
the threshold table reports zero-boundary inference. Neither plot similarity nor
greater precision is interpreted as stronger causal robustness.

### Retained files

The main summary files are:

- `results/main_simulation_summary.csv`
- `results/comparator_summary.csv`
- `results/comparator_seeds.csv`
- `results/comparator_sessionInfo.txt`
- `results/comparator_validation.txt`
- `results/size_power_summary.csv`
- `results/algorithm_sensitivity_summary.csv`
- `results/application_effects.csv`
- `results/application_diagnostics.csv`
- `results/application_covariate_sensitivity.csv`
- `results/application_support.csv`
- `results/application_analysis_specification.txt`
- `results/application_sessionInfo.txt`
- `results/medflex_refresh_sessionInfo.txt`
- `results/paired_fit_summary.csv`
- `results/asset_sessionInfo.txt` (retained-result rebuild only)
- `results/package_versions.csv`
- `results/sessionInfo.txt`
- `results/validation_report.txt`

Numerical failures are retained in the corresponding `*_status.csv` or raw
record file. They are not silently replaced by OLS estimates.
Conditional interval coverage and the fraction of all attempts that both return
an interval and cover the truth are reported separately. Neither summary removes
the need to consider numerical failures when interpreting efficiency gains.

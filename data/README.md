# JOBS II Illustration Data

`jobs.RData` contains the unmodified `jobs` object distributed with version
4.5.1 of the R package `mediation`. It has 899 rows and 17 columns. The object
was checked against the installed package data when this upload snapshot was
prepared. It is not the full original JOBS II study data.

The package documentation describes these data as illustrative, not a basis
for inference about program efficacy. The full original data source is the
ICPSR archive. The manuscript uses this object as a reproducible illustration
and reports model diagnostics and sensitivity analyses.

Analysis variables:

| Role | Variable |
| --- | --- |
| Randomized intervention assignment | `treat` |
| Mediator: job-search self-efficacy | `job_seek` |
| Outcome: post-intervention depressive symptoms | `depress2` |
| Baseline adjustment | `depress1`, `econ_hard`, `sex`, `age` |

Read the data from the repository root with:

```r
load("data/jobs.RData")
```

Upstream documentation is available through `help("jobs", package = "mediation")`.
The package declares `GPL (>= 2)`; these third-party materials are not relicensed
by this repository. The GPL version 2 text is included in
`docs/GPL-2-mediation-data.txt`. Package source: <https://cran.r-project.org/package=mediation>.

References supplied by the upstream package and the manuscript include:

- Vinokur, A. and Schul, Y. (1997). Mastery and inoculation against setbacks as
  active ingredients in the jobs intervention for the unemployed. Journal of
  Consulting and Clinical Psychology, 65(5), 867-877.
- Tingley, D., Yamamoto, T., Hirose, K., Keele, L. and Imai, K. (2014).
  mediation: R Package for Causal Mediation Analysis. Journal of Statistical
  Software, 59(5), 1-38.

The full analysis specification is in `APPLICATION_JOBS2_20260916.md` at the
repository root.

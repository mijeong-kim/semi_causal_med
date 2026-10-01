# Density-score estimation of the Imai sensitivity functional

Target journal: JKSS, as confirmed by the author. No cover-letter changes.
This design is written before fitting the new sensitivity analysis.

## Attribution and model

The sensitivity identification map, rho interpretation, residual-R-squared
calibration and zero-effect boundary are attributed to Imai, Keele and Yamamoto
(2010), DOI 10.1214/10-STS321, Section 5 and Appendix D. This implementation is
not described as the first semiparametric mediation sensitivity analysis.

Fit a pooled mediator regression M on (T,X), retaining its common independent
error distribution, and reduced-form Y regressions on X within each treatment
stratum, not outcome regressions on M. The mediator and reduced-outcome errors
may be dependent; their joint law must not depend on X within a stratum. Their
marginal densities are unspecified. Reduced-outcome coefficients, variances and
cross-covariances may be treatment-specific; mediator variance is common.
A common rho across strata is a sensitivity restriction, not estimated from the
observed data. Conditioning/standardization uses pretreatment X only.

Before the full simulation, the initial four-regression prototype was amended
to pool M. Fully stratifying M would discard the common-shape restriction that
permits density-score efficiency gains for its treatment coefficient. This
amendment follows the influence-function argument, not performance-based design
selection. One Gaussian development data set was used to check the prototype.

For each fixed rho, the causal interpretation additionally requires treatment
exchangeability, positivity, consistency and a linear structural outcome with a
constant mediator slope within each treatment stratum and an additive error
independent of the intervened mediator value. The mediator-outcome no-confounding
assumption is not maintained. This does not accommodate arbitrary nonlinear or
post-treatment confounding.

OLS and the existing marginal density-score regression fits use the same
reduced-form specifications. Stack their scores, the residual cross-product
moments and the empirical baseline mean. Propagate the full joint covariance
through the existing sensitivity map. Marginal efficient scores do not imply
joint semiparametric efficiency under an unrestricted bivariate error law.

The rho=0 density-score reduced-form estimate need not equal the earlier
conditional-outcome density-score estimate in a finite sample. Both outputs must
be labeled; do not silently overwrite or splice the earlier application estimates.

## Fixed numerical evaluation

- n=300, 500 attempted replications per cell, Gaussian/asymmetric-mixture errors.
- True structural correlations: -0.3, 0, 0.3; seed 20260918.
- Main-study structural means and natural effects are unchanged. Generate
  eM=UM and eY=rho*UM+sqrt(1-rho^2)*UY with independent standardized innovations.
- Evaluate pointwise intervals at the true correlation, the zero-correlation
  analysis, and estimates/intervals for the population zero-effect rho threshold.
- Fit three regressions per method (pooled M; reduced Y in each of T=0 and T=1).
- Retain numerical failures and all five effects. Summaries include success,
  bias, RMSE, coverage, MCSE, length, report-and-cover and paired-fit comparisons.
- Additional correlations are alternative scientific assumptions, not values
  chosen to obtain a favorable mediation conclusion.

## Application and plotting

Use the JOBS II primary example and its baseline specification. Produce
PNIE/TNIE curves on rho in [-0.8,0.8], with pointwise nominal 95% intervals.
Produce PNIE residual-variance R-squared contours for BOTH confounding signs,
for OLS and density-score estimates. Axes are residual variance fractions
R_M^{2*}, R_Y^{2*}, NOT the total-variance tilde coordinates in the user's example.
Under Imai's additional variance-decomposition interpretation,
rho=sign*sqrt(R_M^{2*}*R_Y^{2*}). This calibration is not learned from the data.
Distinguish the point-estimate-zero contour from pointwise interval-zero contours.
No claim of simultaneous confidence bands, verified ignorability, general
heteroscedastic robustness, or empirical superiority is made.

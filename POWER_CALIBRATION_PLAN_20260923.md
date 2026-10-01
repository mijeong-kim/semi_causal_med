# Independent null calibration and Huber power comparison

Specification written before generating the new Monte Carlo results.
This supplements, and does not replace, the original size-power experiment.

- Target: PNIE = beta2 * gamma; two-sided testing of PNIE = 0.
- Same original asymmetric-mixture DGP, n=300, gamma=-0.26, eta=0.8;
  all other coefficients and the error distribution unchanged.
- Methods: existing OLS and semiparametric estimators, Huber-FIX and Huber-SEL.
- Calibration: 5,000 attempted independent null data sets, beta2=0.
- Evaluation: 1,000 new data sets at each beta2 in 0, 0.1, 0.2, 0.3, 0.4.
- Master seed 20260923 generates unique task seeds for all 10,000 data sets.
  Every method analyzes the same data within a task; no evaluation data enter
  threshold calibration. The held-out null is not the calibration sample.
- Critical value: the ceil(0.95*(m+1))-th ordered absolute Wald statistic
  among m successful null-calibration fits. Reject strictly above it.
  This is a DGP-specific diagnostic, not a general-purpose calibrated test
  or a uniform-size guarantee over the composite mediation null.
- Report nominal and calibrated rejection, valid-fit counts, MCSE, numerical
  success and all-attempt report-and-reject fractions. Also report calibrated
  contrasts on exactly the same successful evaluation data sets.
- Resample calibration tasks and evaluation tasks independently, pairing methods
  within each resample, for 1,000 bootstrap summaries (seed 20260924). Report
  uncertainty from both simulation stages, not just conditional binomial MCSE.
- Preserve all attempts, seeds, failure messages, tuning parameters and source
  hashes. Source-checked checkpoints permit resumption without changing seeds.

## Huber comparator specification

These are interaction-adapted Huber slope-product comparators motivated by
Wang, Peng and Tong (2025), DOI 10.1017/psy.2024.28. They are not a verbatim
replication of that paper's simple-mediation Sobel implementation.

Both regressions include an intercept and the same terms as the other methods:
M ~ T + X and Y ~ T*M + X. Obtain a LAD pilot using quantreg::rq, tau=0.5.
Let s = median(abs(LAD residual))/0.6745. Huber-FIX uses raw threshold 1.345*s;
Huber-SEL minimizes mean(min(e^2,k^2))/mean(abs(e)<=k)^2 on the author's grid
0.2,0.21,...,round(3*s,2), using the LAD residuals. If the upper endpoint is
less than 0.2, report failure rather than silently changing the grid.
Estimate coefficients by convex Huber IRLS, with the threshold held fixed.

Use final-residual, centered-design slope influence functions and their joint
empirical covariance to propagate the PNIE product. Under the independent,
common error law, the intercept absorbs the constant Huber-location shift;
slopes and therefore PNIE target the same parameters even for asymmetric errors.
The centered-design slope equations are first-order orthogonal to the intercept
and threshold, under that same independence restriction. Do not use Huber
intercepts as mean intercepts for direct or total effects in this experiment.
This covariance adaptation retains cross-regression dependence and evaluates
Huber scores at the fitted location instead of the LAD pilot location.

The authors' public Ck.R and Inference.R were inspected on 2026-09-23 at
https://github.com/pxj66/REMA (GPL-3.0). The implementation here is written
from the mathematical specification; no author code is vendored.

# Full-curve extension of the retained confounding simulation

The author requested that the title reflect semiparametric inference and
sensitivity analysis and that sensitivity analysis also be performed in the
simulation. This extension is specified before examining full-curve results.
It does not select new favorable data sets or alter the estimators.

## Design

- Refit the SAME 3,000 data sets in the existing confounding experiment: n=300,
  500 replications, Gaussian/asymmetric-mixture innovations and true structural
  correlations -0.3, 0, 0.3. Recreate the original task grid with seed 20260918.
- Compare the same OLS-RF and DS-RF estimators. Do not change the density
  bandwidth, starts, root-selection rule, working models or failure criteria.
- Evaluate the common assumed correlation on 33 points from -0.8 to 0.8 in
  increments of 0.05. All five effects are retained, with PNIE/TNIE emphasized
  in figures. These are paired evaluations, not new independent replications.
- Retain failed attempts and messages. Verify the repeated fits against the
  original estimates, standard errors and zero boundaries before interpretation.

## Two distinct targets

For a true generating correlation r and an assumed value rho, the population
sensitivity slope is b_t+r-h(rho)*sqrt(1-r^2), where b_t=-0.8+t and
h(rho)=rho/sqrt(1-rho^2). Its indirect effect is 0.4 times that slope.
The population total effect is 0.78; obtain direct effects by subtraction.

Curve bias, RMSE and nominal pointwise coverage target this population
sensitivity functional at the specified rho. Separately retain bias and
interval inclusion for the generating causal effects (-0.32,0.08,0.70,1.10,0.78).
The targets coincide at rho=r; away from r, causal-effect inclusion is not a
coverage-calibration criterion for the assumed-rho functional.

## Outputs and interpretation

Under a common rho each estimated effect has the exact form A+h(rho)*B, with
variance V_A+2*h(rho)*C_AB+h(rho)^2*V_B. Store these five coefficients for every
effect and attempted fit, plus status, seeds and refitted zero boundaries.
This compact representation reconstructs every pointwise estimate and interval
without storing nearly one million redundant grid rows or refitting again.

Provide complete grid summaries with MCSEs, all-attempted report-and-cover,
paired-fit ratios and paired coverage. Main-article plots show population curves,
Monte Carlo mean estimates and central 90% empirical ranges across successful
replications. The latter describe sampling variation, not confidence bands.
Supplementary plots show curve-functional coverage and paired RMSE ratios for
PNIE/TNIE across every generating correlation and both innovation densities.
Summarize zero-boundary bias, RMSE, coverage and interval length for both effects
from the retained refits; no improvement is assumed in advance.

No new causal identification formula, simultaneous band, joint efficiency bound,
or correction for variance-model misspecification is claimed. Preserve all
previous numerical files, the original primary application and the cover letter.

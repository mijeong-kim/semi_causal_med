# Semiparametric Causal Mediation Analysis: R Code

This repository is intended as a code release for the semiparametric causal mediation project.

GitHub repository: https://github.com/mijeong-kim/semi_causal_med

The public GitHub version contains R code only. Manuscript source files and compiled PDFs are not included in the code repository.

## Included R scripts

- `scripts/jci_interaction_mediation_revision.R`
  - simulation study
  - near-boundary power design
  - `jobs` application
  - `uis` application
- `scripts/jci_interaction_figures.R`
  - figure generation for the empirical applications

## Suggested GitHub contents

For a code-only public release, upload at least:

- `scripts/`
- `README.md`

Optional additions:

- `results/` if you want to share derived numerical outputs
- `figures/` if you want to share manuscript figure files generated from the code

## Notes

- The `uis` analysis uses the rescaled outcome `TIME/100` only for numerical conditioning in root finding.
- The code corresponds to the semiparametric mediation manuscript submitted to the *Journal of Causal Inference*.

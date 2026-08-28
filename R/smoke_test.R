source(file.path("R", "semiparametric_mediation.R"))

required_cs_packages(include_comparators = TRUE)
set.seed(20260827)
data <- generate_mediation_data(120, "asymmetric_mixture")

ols <- fit_stacked_mediation(data, "OLS", baseline_names = "X")
proposed <- fit_stacked_mediation(
  data,
  "Semiparametric",
  baseline_names = "X"
)
quasi <- fit_mediation_package(
  data,
  boot = FALSE,
  sims = 50,
  seed = 20260828
)
bootstrap <- suppressMessages(fit_mediation_package(
  data,
  boot = TRUE,
  sims = 49,
  seed = 20260829
))
nem <- suppressMessages(fit_medflex(data))

expected_effects <- c("PNIE", "TNIE", "PNDE", "TNDE", "TE")
outputs <- list(ols, proposed, quasi, bootstrap, nem)
stopifnot(all(vapply(outputs, function(x) {
  identical(as.character(x$Effect), expected_effects) &&
    all(is.finite(x$Estimate)) &&
    all(is.finite(x$Lower)) &&
    all(is.finite(x$Upper))
}, logical(1))))

decomposition_error <- vapply(outputs, function(x) {
  effect <- stats::setNames(x$Estimate, x$Effect)
  max(abs(c(
    effect["TE"] - effect["PNDE"] - effect["TNIE"],
    effect["TE"] - effect["TNDE"] - effect["PNIE"]
  )))
}, numeric(1))
stopifnot(all(decomposition_error < 1e-8))

cat(
  "Smoke test passed for all estimators and both total-effect decompositions.\n"
)

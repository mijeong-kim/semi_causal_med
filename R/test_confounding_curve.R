source(file.path("R", "confounding_curve.R"))
truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.70, TNDE = 1.10, TE = 0.78)
for (r in c(-0.3, 0, 0.3)) {
  theta <- c("m:(Intercept)" = 0.2, "m:T" = 0.4, "m:X" = 0.3, "m:log_sigma2" = 0,
             "w0:(Intercept)" = -0.16, "w0:X" = 0.16,
             "w0:log_sigma2" = log((-0.8 + r)^2 + 1 - r^2),
             "w1:(Intercept)" = 0.62, "w1:X" = 0.46,
             "w1:log_sigma2" = log((0.2 + r)^2 + 1 - r^2),
             cross0 = -0.8 + r, cross1 = 0.2 + r, "xbar:X" = 0)
  for (rho in curve_rho_grid()) {
    stopifnot(max(abs(curve_population(r, rho) - sensitivity_map(theta, "X", rho)$estimate)) < 1e-12)
  }
  stopifnot(max(abs(curve_population(r, r) - truth)) < 1e-12)
}
set.seed(20260919)
d <- confounding_data(300, "gaussian", 0.3)
max_error <- 0
for (method in c("OLS", "Semiparametric")) {
  fit <- fit_reduced_sensitivity(d, method, "X")
  basis <- curve_basis(fit)
  stopifnot(all(basis$VA > 0), all(basis$VB >= 0),
            all(basis$CAB^2 <= basis$VA * basis$VB + 1e-12),
            abs(basis$B[basis$Effect == "TE"]) < 1e-12)
  for (rho in curve_rho_grid()) {
    direct <- sensitivity_effects(fit, rho)
    compact <- curve_evaluate(basis, rho)
    max_error <- max(max_error, abs(as.matrix(compact) - as.matrix(direct[names(compact)])))
    b <- setNames(compact$Estimate, basis$Effect)
    stopifnot(abs(b["TE"] - b["PNIE"] - b["TNDE"]) < 1e-12,
              abs(b["TE"] - b["TNIE"] - b["PNDE"]) < 1e-12)
  }
}
stopifnot(max_error < 1e-10)
cat("Full-curve population, compact-covariance and decomposition tests PASSED. Maximum error:",
    format(max_error), "\n")

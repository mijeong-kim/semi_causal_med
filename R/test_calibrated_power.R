source(file.path("R", "calibrated_power.R"))
set.seed(9301)
data <- generate_mediation_data(300L, "asymmetric_mixture", beta2 = 0.1,
                                gamma = -0.26, eta = 0.8)
ols <- fit_ols_lm(M ~ T + X, data)
large <- fit_huber_slopes(M ~ T + X, data, threshold = 1e6)
stopifnot(max(abs(ols$coefficients - large$coefficients)) < 1e-8,
          max(abs(crossprod(large$influence) / nrow(data)^2 - ols$covariance[-1, -1])) < 1e-8)
large_out <- fit_huber_slopes(Y ~ T * M + X, data, threshold = 1e6)
joint_if <- large_out$coefficients["M"] * large$influence[, "T"] +
  large$coefficients["T"] * large_out$influence[, "M"]
ols_effects <- fit_stacked_mediation(data, "OLS", "X")
stopifnot(abs(sqrt(sum(joint_if^2)) / nrow(data) -
              ols_effects$StdError[ols_effects$Effect == "PNIE"]) < 1e-8)
for (selected in c(FALSE, TRUE)) {
  fit <- fit_huber_slopes(Y ~ T * M + X, data, selected)
  shifted <- data
  shifted$Y <- shifted$Y + 17
  other <- fit_huber_slopes(Y ~ T * M + X, shifted, selected)
  stopifnot(max(abs(fit$coefficients[-1] - other$coefficients[-1])) < 1e-7,
            abs(other$coefficients[1] - fit$coefficients[1] - 17) < 1e-7,
            max(abs(fit$influence - other$influence)) < 1e-7)
  # Independent optimizer checks the IRLS solution of the convex Huber loss.
  design <- model.matrix(Y ~ T * M + X, data)
  threshold <- fit$threshold
  objective <- function(beta) {
    error <- abs(data$Y - design %*% beta)
    sum(ifelse(error <= threshold, error^2 / 2, threshold * error - threshold^2 / 2))
  }
  reference <- optim(coef(lm(Y ~ T * M + X, data)), objective, method = "BFGS",
                     control = list(reltol = 1e-12, maxit = 2000))
  stopifnot(reference$convergence == 0,
            abs(reference$value - objective(fit$coefficients)) < 1e-7,
            max(abs(reference$par - fit$coefficients)) < 1e-4)
  pnie <- fit_huber_pnie(data, selected)
  stopifnot(is.finite(pnie$Estimate), pnie$StdError > 0)
}
stopifnot(calibrated_critical(1:100) == 96,
          calibrated_critical(c(1:100, NA, Inf)) == 96)
grid <- power_seed_grid()
stopifnot(nrow(grid) == 10000, !anyDuplicated(grid$Seed),
          sum(grid$Stage == "calibration") == 5000,
          all(table(grid$Beta2[grid$Stage == "evaluation"]) == 1000))
message("Huber coefficients, covariance, shift invariance, optimizer and calibration tests passed.")

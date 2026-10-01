source(file.path("R", "confounding_sensitivity.R"))
set.seed(88216)
data <- confounding_data(600, "gaussian", 0.3)
fit <- fit_reduced_sensitivity(data, "OLS", "X")
theta <- fit$parameter
max_error <- 0
for (rho in c(-0.7, -0.3, 0, 0.3, 0.7)) {
  mapped <- sensitivity_map(theta, "X", rho)
  numeric <- sapply(seq_along(theta), function(j) {
    h <- 1e-6 * max(1, abs(theta[j]))
    increment <- numeric(length(theta))
    increment[j] <- h
    (sensitivity_map(theta + increment, "X", rho)$estimate -
       sensitivity_map(theta - increment, "X", rho)$estimate) / (2 * h)
  })
  max_error <- max(max_error, abs(mapped$gradient - numeric))
  stopifnot(abs(mapped$estimate["TE"] - mapped$estimate["PNDE"] - mapped$estimate["TNIE"]) < 1e-12,
            abs(mapped$estimate["TE"] - mapped$estimate["TNDE"] - mapped$estimate["PNIE"]) < 1e-12)
}
stopifnot(max_error < 1e-7, max(abs(colMeans(fit$score))) < 1e-10,
          max(abs(fit$covariance - t(fit$covariance))) < 1e-12,
          min(eigen(fit$covariance, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)

# Reproduce Imai's covariance identity independently, including at population truth.
for (rho_true in c(-0.3, 0, 0.3)) for (t in 0:1) {
  b <- -0.8 + t
  vm <- 1
  vw <- (b + rho_true)^2 + 1 - rho_true^2
  cv <- b + rho_true
  r <- cv / sqrt(vm * vw)
  b_recovered <- sqrt(vw / vm) * (r - rho_true * sqrt((1 - r^2) / (1 - rho_true^2)))
  stopifnot(abs(b_recovered - b) < 1e-12)
}
for (t in 0:1) {
  z <- sensitivity_map(theta, "X", 0)
  r <- z$rho_zero[t + 1]
  value <- sensitivity_map(theta, "X", r)$estimate[c("PNIE", "TNIE")[t + 1]]
  stopifnot(abs(value) < 1e-12)
  numerical <- sapply(seq_along(theta), function(j) {
    h <- 1e-6
    inc <- numeric(length(theta)); inc[j] <- h
    (sensitivity_map(theta + inc, "X", 0)$rho_zero[t + 1] -
       sensitivity_map(theta - inc, "X", 0)$rho_zero[t + 1]) / (2 * h)
  })
  stopifnot(max(abs(numerical - z$rho_gradient[t + 1, ])) < 1e-7)
  rows <- data$T == t
  em <- resid(lm(M ~ T + X, data))[rows]
  ew <- resid(lm(Y ~ X, data[rows, ]))
  independent <- unname(coef(lm(M ~ T + X, data))["T"] * mean(em * ew) /
                         mean(resid(lm(M ~ T + X, data))^2))
  stopifnot(abs(independent - z$estimate[t + 1]) < 1e-12)
}

# The same covariance-to-slope map also works with separate rho_0 and rho_1.
pair <- sensitivity_map(theta, "X", c(-0.2, 0.4))
stopifnot(abs(pair$estimate["PNIE"] - sensitivity_map(theta, "X", -0.2)$estimate["PNIE"]) < 1e-12,
          abs(pair$estimate["TNIE"] - sensitivity_map(theta, "X", 0.4)$estimate["TNIE"]) < 1e-12)

cat("Confounding sensitivity algebra/gradient/stack tests PASSED. Maximum gradient error:", max_error, "\n")

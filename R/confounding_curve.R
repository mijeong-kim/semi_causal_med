source(file.path("R", "confounding_sensitivity.R"))

curve_rho_grid <- function() seq(-16, 16) / 20

curve_population <- function(true_rho, assumed_rho) {
  stopifnot(abs(true_rho) < 1, abs(assumed_rho) < 1)
  h <- assumed_rho / sqrt(1 - assumed_rho^2)
  delta <- 0.4 * (c(-0.8, 0.2) + true_rho - h * sqrt(1 - true_rho^2))
  c(PNIE = delta[1], TNIE = delta[2], PNDE = 0.78 - delta[2],
    TNDE = 0.78 - delta[1], TE = 0.78)
}

curve_basis <- function(fit) {
  zero <- sensitivity_map(fit$parameter, fit$baseline_names, 0)
  unit <- sensitivity_map(fit$parameter, fit$baseline_names, 1 / sqrt(2))
  ga <- zero$gradient
  gb <- unit$gradient - ga
  v <- fit$covariance
  data.frame(Effect = names(zero$estimate), A = unname(zero$estimate),
             B = unname(unit$estimate - zero$estimate),
             VA = unname(diag(ga %*% v %*% t(ga))),
             CAB = unname(diag(ga %*% v %*% t(gb))),
             VB = unname(diag(gb %*% v %*% t(gb))))
}

curve_evaluate <- function(basis, rho) {
  stopifnot(length(rho) == 1L, is.finite(rho), abs(rho) < 1)
  h <- rho / sqrt(1 - rho^2)
  estimate <- basis$A + h * basis$B
  variance <- basis$VA + 2 * h * basis$CAB + h^2 * basis$VB
  if (any(variance <= 0, na.rm = TRUE)) stop("Nonpositive curve variance.")
  se <- sqrt(variance)
  data.frame(Estimate = estimate, StdError = se,
             Lower = estimate - qnorm(0.975) * se,
             Upper = estimate + qnorm(0.975) * se)
}

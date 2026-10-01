# PNIE-only Huber comparison and independent-null calibration utilities.
source(file.path("R", "semiparametric_mediation.R"))

power_methods <- c("OLS", "Huber-FIX", "Huber-SEL", "Semiparametric")

fit_huber_slopes <- function(formula, data, selected = FALSE, threshold = NULL) {
  if (!requireNamespace("quantreg", quietly = TRUE)) stop("Install quantreg.")
  frame <- stats::model.frame(formula, data)
  design <- stats::model.matrix(formula, frame)
  response <- stats::model.response(frame)
  stopifnot(colnames(design)[1L] == "(Intercept)", all(design[, 1L] == 1))
  pilot <- quantreg::rq(formula, data = data, tau = 0.5, method = "br")
  residual <- as.numeric(stats::residuals(pilot))
  scale <- stats::median(abs(residual)) / 0.6745
  if (!is.finite(scale) || scale <= 0) stop("Invalid Huber pilot scale.")
  if (is.null(threshold)) {
    threshold <- 1.345 * scale
    if (selected) {
      upper <- round(3 * scale, 2)
      if (upper < 0.2) stop("Empty Huber-SEL tuning grid.")
      grid <- seq(0.2, upper, by = 0.01)
      tau <- vapply(grid, function(k) {
        fraction <- mean(abs(residual) <= k)
        if (fraction == 0) return(Inf)
        mean(pmin(residual^2, k^2)) / fraction^2
      }, numeric(1))
      threshold <- grid[which.min(tau)]
    }
  }
  stopifnot(is.finite(threshold), threshold > 0)
  coefficient <- stats::coef(pilot)
  converged <- FALSE
  for (iteration in seq_len(500L)) {
    residual <- as.numeric(response - design %*% coefficient)
    weight <- pmin(1, threshold / pmax(abs(residual), .Machine$double.eps))
    updated <- stats::lm.wfit(design, response, weight)$coefficients
    if (any(!is.finite(updated))) stop("Singular Huber weighted design.")
    converged <- max(abs(updated - coefficient)) < 1e-10 * (1 + max(abs(coefficient)))
    coefficient <- updated
    if (converged) break
  }
  if (!converged) stop("Huber IRLS did not converge.")
  residual <- as.numeric(response - design %*% coefficient)
  psi <- pmax(-threshold, pmin(threshold, residual))
  if (max(abs(colMeans(design * psi))) > 1e-7) stop("Huber score tolerance failed.")
  centered <- scale(design[, -1L, drop = FALSE], center = TRUE, scale = FALSE)
  derivative <- mean(abs(residual) < threshold)
  gram <- crossprod(centered) / nrow(design)
  # Independence makes slope equations orthogonal to the intercept and threshold.
  influence <- (centered %*% safe_inverse(gram)) * psi / derivative
  influence <- scale(influence, center = TRUE, scale = FALSE)
  if (any(!is.finite(influence))) stop("Invalid Huber slope influence.")
  list(coefficients = coefficient, influence = influence, threshold = threshold,
       scale = scale, iterations = iteration)
}

fit_huber_pnie <- function(data, selected = FALSE) {
  med <- fit_huber_slopes(M ~ T + X, data, selected)
  out <- fit_huber_slopes(Y ~ T * M + X, data, selected)
  a <- unname(med$coefficients["T"])
  b <- unname(out$coefficients["M"])
  influence <- b * med$influence[, "T"] + a * out$influence[, "M"]
  list(Estimate = a * b, StdError = sqrt(sum(influence^2)) / nrow(data),
       ThresholdM = med$threshold, ThresholdY = out$threshold)
}

power_seed_grid <- function(calibration = 5000L, evaluation = 1000L) {
  evaluation_grid <- expand.grid(Replication = seq_len(evaluation),
                                 Beta2 = c(0, 0.1, 0.2, 0.3, 0.4))
  grid <- rbind(data.frame(Stage = "calibration", Replication = seq_len(calibration), Beta2 = 0),
                data.frame(Stage = "evaluation", evaluation_grid))
  set.seed(20260923)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid), replace = FALSE)
  grid$Task <- seq_len(nrow(grid))
  grid
}

calibrated_critical <- function(statistic, alpha = 0.05) {
  statistic <- statistic[is.finite(statistic)]
  if (!length(statistic)) stop("No valid calibration statistics.")
  rank <- ceiling((1 - alpha) * (length(statistic) + 1))
  if (rank > length(statistic)) return(Inf)
  sort(statistic, partial = rank)[rank]
}

run_power_task <- function(task) {
  set.seed(task$Seed)
  data <- generate_mediation_data(300L, "asymmetric_mixture", beta2 = task$Beta2,
                                  gamma = -0.26, eta = 0.8)
  do.call(rbind, lapply(power_methods, function(method) {
    warnings <- character()
    error <- ""
    fit <- tryCatch(withCallingHandlers({
      if (startsWith(method, "Huber")) {
        fit_huber_pnie(data, selected = method == "Huber-SEL")
      } else {
        full <- fit_stacked_mediation(data, method, "X")
        as.list(full[full$Effect == "PNIE", c("Estimate", "StdError")])
      }
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) {
      error <<- conditionMessage(e)
      NULL
    })
    valid <- !is.null(fit) && is.finite(fit$Estimate) &&
      is.finite(fit$StdError) && fit$StdError > 0
    data.frame(task, Method = method, Success = valid,
      TruePNIE = -0.26 * task$Beta2,
      Estimate = if (valid) fit$Estimate else NA_real_,
      StdError = if (valid) fit$StdError else NA_real_,
      AbsZ = if (valid) abs(fit$Estimate / fit$StdError) else NA_real_,
      ThresholdM = if (valid && !is.null(fit$ThresholdM)) fit$ThresholdM else NA_real_,
      ThresholdY = if (valid && !is.null(fit$ThresholdY)) fit$ThresholdY else NA_real_,
      Error = error, Warning = paste(unique(warnings), collapse = " | "),
      stringsAsFactors = FALSE)
  }))
}

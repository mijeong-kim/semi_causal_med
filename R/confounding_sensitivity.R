source(file.path("R", "semiparametric_mediation.R"))

fit_ols_reduced <- function(formula, data) {
  fit <- fit_ols_lm(formula, data)
  residual <- drop(fit$response - fit$design %*% fit$coefficients)
  v <- mean(residual^2)
  p <- length(fit$par)
  fit$par <- c(fit$par, log_sigma2 = log(v))
  fit$score <- cbind(fit$score, log_sigma2 = residual^2 / v - 1)
  bread <- matrix(0, p + 1, p + 1)
  bread[seq_len(p), seq_len(p)] <- fit$bread
  bread[p + 1, seq_len(p)] <- -2 * colMeans(fit$design * residual) / v
  bread[p + 1, p + 1] <- -1
  fit$bread <- bread
  fit
}

fit_reduced_sensitivity <- function(data, method = c("Semiparametric", "OLS"),
                                    baseline_names = character()) {
  method <- match.arg(method)
  stopifnot(all(data$T %in% 0:1), !anyNA(data[c("M", "Y", "T", baseline_names)]))
  fitter <- if (method == "OLS") fit_ols_reduced else fit_semiparametric_lm
  fits <- list(m = fitter(reformulate(c("T", baseline_names), "M"), data))
  for (t in 0:1) {
    formula <- if (length(baseline_names)) reformulate(baseline_names, "Y") else Y ~ 1
    fits[[paste0("w", t)]] <- fitter(formula, data[data$T == t, , drop = FALSE])
  }
  n <- nrow(data)
  parameter <- unlist(lapply(names(fits), function(key) {
    setNames(fits[[key]]$par, paste0(key, ":", names(fits[[key]]$par)))
  }))
  crosses <- sapply(0:1, function(t) {
    m <- fits$m
    w <- fits[[paste0("w", t)]]
    mean((m$response - drop(m$design %*% m$coefficients))[data$T == t] *
           (w$response - drop(w$design %*% w$coefficients)))
  })
  parameter <- c(parameter, setNames(crosses, c("cross0", "cross1")))
  if (length(baseline_names)) {
    parameter <- c(parameter, setNames(colMeans(data[baseline_names]), paste0("xbar:", baseline_names)))
  }
  p <- length(parameter)
  score <- matrix(0, n, p, dimnames = list(NULL, names(parameter)))
  bread <- matrix(0, p, p, dimnames = list(names(parameter), names(parameter)))
  columns <- paste0("m:", names(fits$m$par))
  score[, columns] <- fits$m$score
  bread[columns, columns] <- fits$m$bread
  for (t in 0:1) {
    rows <- which(data$T == t)
    key <- paste0("w", t)
    fit <- fits[[key]]
    columns <- paste0(key, ":", names(fit$par))
    score[rows, columns] <- fit$score
    bread[columns, columns] <- length(rows) / n * fit$bread
    m <- fits$m
    w <- fits[[paste0("w", t)]]
    em <- (m$response - drop(m$design %*% m$coefficients))[rows]
    ew <- w$response - drop(w$design %*% w$coefficients)
    cross_name <- paste0("cross", t)
    score[rows, cross_name] <- em * ew - parameter[cross_name]
    bread[cross_name, cross_name] <- -length(rows) / n
    bread[cross_name, paste0("m:", names(m$coefficients))] <- -colSums(m$design[rows, , drop = FALSE] * ew) / n
    bread[cross_name, paste0("w", t, ":", names(w$coefficients))] <- -colSums(w$design * em) / n
  }
  if (length(baseline_names)) {
    columns <- paste0("xbar:", baseline_names)
    score[, columns] <- sweep(as.matrix(data[baseline_names]), 2, parameter[columns], "-")
    bread[columns, columns] <- -diag(length(columns))
  }
  inverse <- safe_inverse(bread)
  covariance <- inverse %*% crossprod(score) %*% t(inverse) / n^2
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) stop("Invalid reduced-form covariance.")
  result <- list(parameter = parameter, covariance = covariance, score = score, bread = bread,
                 fits = fits, baseline_names = baseline_names, method = method)
  sensitivity_map(parameter, baseline_names, 0)
  result
}

sensitivity_map <- function(theta, baseline_names = character(), rho = 0) {
  if (length(rho) == 1L) rho <- rep(rho, 2)
  stopifnot(length(rho) == 2, all(is.finite(rho)), all(abs(rho) < 1))
  p <- length(theta)
  zeros <- function() setNames(numeric(p), names(theta))
  reference <- c("(Intercept)" = 1, setNames(theta[paste0("xbar:", baseline_names)], baseline_names))
  contrast <- function(prefix) {
    coefficients <- theta[paste0(prefix, "1:", names(reference))] -
      theta[paste0(prefix, "0:", names(reference))]
    value <- sum(coefficients * reference)
    gradient <- zeros()
    gradient[paste0(prefix, "1:", names(reference))] <- reference
    gradient[paste0(prefix, "0:", names(reference))] <- -reference
    if (length(baseline_names)) gradient[paste0("xbar:", baseline_names)] <- coefficients[-1]
    list(value = value, gradient = gradient)
  }
  dm <- list(value = unname(theta["m:T"]), gradient = zeros())
  dm$gradient["m:T"] <- 1
  total <- contrast("w")
  delta <- numeric(2)
  gradients <- matrix(0, 2, p)
  correlations <- numeric(2)
  correlation_gradients <- matrix(0, 2, p)
  for (t in 0:1) {
    vm_name <- "m:log_sigma2"
    vw_name <- paste0("w", t, ":log_sigma2")
    c_name <- paste0("cross", t)
    vm <- exp(theta[vm_name])
    vw <- exp(theta[vw_name])
    cv <- theta[c_name]
    slope0 <- cv / vm
    d2 <- vw / vm - slope0^2
    if (!is.finite(d2) || d2 <= 0) stop("Nonpositive reduced-form conditional residual variance.")
    d <- sqrt(d2)
    h <- rho[t + 1] / sqrt(1 - rho[t + 1]^2)
    slope <- slope0 - h * d
    db <- zeros()
    db[c_name] <- (1 + h * slope0 / d) / vm
    db[vm_name] <- -slope0 + h * (vw / vm - 2 * slope0^2) / (2 * d)
    db[vw_name] <- -h * vw / vm / (2 * d)
    delta[t + 1] <- dm$value * slope
    gradients[t + 1, ] <- slope * dm$gradient + dm$value * db
    correlations[t + 1] <- cv / sqrt(vm * vw)
    dr <- zeros()
    dr[c_name] <- 1 / sqrt(vm * vw)
    dr[c(vm_name, vw_name)] <- -correlations[t + 1] / 2
    correlation_gradients[t + 1, ] <- dr
  }
  estimate <- c(PNIE = delta[1], TNIE = delta[2], PNDE = total$value - delta[2],
                TNDE = total$value - delta[1], TE = total$value)
  gradient <- rbind(gradients, total$gradient - gradients[2, ],
                    total$gradient - gradients[1, ], total$gradient)
  dimnames(gradient) <- list(names(estimate), names(theta))
  dimnames(correlation_gradients) <- list(c("PNIE", "TNIE"), names(theta))
  list(estimate = estimate, gradient = gradient, mean_m_difference = dm$value,
       rho_zero = correlations, rho_gradient = correlation_gradients)
}

sensitivity_effects <- function(fit, rho = 0) {
  values <- sensitivity_map(fit$parameter, fit$baseline_names, rho)
  se <- sqrt(diag(values$gradient %*% fit$covariance %*% t(values$gradient)))
  if (any(!is.finite(se)) || any(se <= 0)) stop("Invalid sensitivity standard error.")
  data.frame(Effect = names(values$estimate), Estimate = unname(values$estimate),
             StdError = unname(se), Lower = unname(values$estimate - qnorm(0.975) * se),
             Upper = unname(values$estimate + qnorm(0.975) * se), Method = fit$method)
}

sensitivity_thresholds <- function(fit) {
  values <- sensitivity_map(fit$parameter, fit$baseline_names, 0)
  r <- values$rho_zero
  se <- sqrt(diag(values$rho_gradient %*% fit$covariance %*% t(values$rho_gradient)))
  zse <- se / (1 - r^2)
  lower <- tanh(atanh(r) - qnorm(0.975) * zse)
  upper <- tanh(atanh(r) + qnorm(0.975) * zse)
  data.frame(Effect = c("PNIE", "TNIE"), RhoZero = r, StdError = se,
             Lower = lower, Upper = upper, ResidualR2Product = r^2,
             MeanMediatorDifference = values$mean_m_difference, Method = fit$method)
}

confounding_data <- function(n, scenario, rho) {
  X <- rnorm(n)
  T <- rbinom(n, 1, 0.5)
  um <- generate_standardized_error(n, scenario)
  uy <- generate_standardized_error(n, scenario)
  M <- 0.2 + 0.4 * T + 0.3 * X + um
  Y <- 0.5 * T - 0.8 * M + T * M + 0.4 * X + rho * um + sqrt(1 - rho^2) * uy
  data.frame(Y = Y, M = M, T = T, X = X)
}

confounding_truth <- function(rho) {
  slopes <- c(-0.8, 0.2)
  covariance <- slopes + rho
  covariance / sqrt(covariance^2 + 1 - rho^2)
}

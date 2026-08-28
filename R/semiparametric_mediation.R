required_cs_packages <- function(include_comparators = FALSE) {
  packages <- c("Matrix", "rootSolve")
  if (include_comparators) {
    packages <- c(packages, "mediation", "medflex")
  }

  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

safe_inverse <- function(mat, rcond_tol = 1e-10) {
  rc <- suppressWarnings(rcond(mat))
  if (!is.finite(rc) || rc < rcond_tol) {
    stop("Matrix is singular or nearly singular (rcond < ", rcond_tol, ").")
  }
  solve(mat)
}

semi_score_i <- function(par, design, response, density_floor = 1e-10) {
  p <- ncol(design)
  beta <- par[seq_len(p)]
  sigma2 <- exp(par[p + 1L])
  residual <- drop(response - design %*% beta)
  n <- length(residual)

  bandwidth <- 1.06 * stats::sd(residual) * n^(-1 / 5)
  if (!is.finite(bandwidth) || bandwidth <= 0) {
    stop("The residual bandwidth is not positive and finite.")
  }

  differences <- outer(residual, residual, "-")
  kernel <- stats::dnorm(differences / bandwidth)
  diag(kernel) <- 0
  density <- rowSums(kernel) / ((n - 1) * bandwidth)
  density_derivative <- rowSums(-differences * kernel) /
    ((n - 1) * bandwidth^3)
  density_score <- density_derivative / pmax(density, density_floor)

  mu3 <- mean(residual^3)
  mu4 <- mean(residual^4)
  q_residual <- residual^2 - sigma2 - mu3 * residual / sigma2
  eq2 <- mu4 - sigma2^2 - mu3^2 / sigma2
  if (!is.finite(eq2) || eq2 <= density_floor) {
    stop("The estimated second moment of q(epsilon) is not positive.")
  }

  design_mean <- colMeans(design)
  centered_design <- sweep(design, 2, design_mean, "-")
  location_term <- residual / sigma2 -
    mu3 * q_residual / (sigma2 * eq2)

  score_beta <- -centered_design * density_score +
    outer(location_term, design_mean)
  score_log_sigma2 <- sigma2 * q_residual / eq2

  score <- cbind(score_beta, score_log_sigma2)
  colnames(score) <- c(colnames(design), "log_sigma2")
  score
}

semi_score_mean <- function(par, design, response) {
  tryCatch(
    colMeans(semi_score_i(par, design, response)),
    error = function(e) rep(1e6, length(par))
  )
}

numeric_bread <- function(par, score_function, rel_step = 1e-4) {
  p <- length(par)
  bread <- matrix(NA_real_, p, p)
  for (j in seq_len(p)) {
    step_size <- rel_step * max(1, abs(par[j]))
    step <- rep(0, p)
    step[j] <- step_size
    bread[, j] <- (score_function(par + step) - score_function(par - step)) /
      (2 * step_size)
  }
  bread
}

make_start_candidates <- function(ols_beta, ols_sigma2, step_scale = 0.10) {
  coefficient_step <- pmax(0.05, step_scale * pmax(1, abs(ols_beta)))
  alternating <- rep(c(1, -1), length.out = length(ols_beta))
  base <- c(ols_beta, log(ols_sigma2))

  list(
    base,
    c(ols_beta + coefficient_step, log(ols_sigma2)),
    c(ols_beta - coefficient_step, log(ols_sigma2)),
    c(ols_beta + alternating * coefficient_step, log(ols_sigma2)),
    c(ols_beta - alternating * coefficient_step, log(ols_sigma2)),
    c(ols_beta, log(1.25 * ols_sigma2)),
    c(ols_beta, log(0.80 * ols_sigma2))
  )
}

root_distance_from_ols <- function(par, ols_beta, ols_sigma2, design) {
  p <- length(ols_beta)
  beta_difference <- par[seq_len(p)] - ols_beta
  gram <- crossprod(design) / nrow(design)
  coefficient_distance <- drop(
    t(beta_difference) %*% gram %*% beta_difference / ols_sigma2
  )
  variance_distance <- (par[p + 1L] - log(ols_sigma2))^2
  coefficient_distance + variance_distance
}

fit_semiparametric_lm <- function(
  formula,
  data,
  maxiter = 200,
  step_scale = 0.10,
  score_tol = 1e-6,
  rcond_tol = 1e-10
) {
  required_cs_packages()
  ols <- stats::lm(formula, data = data)
  model_frame <- stats::model.frame(ols)
  response <- stats::model.response(model_frame)
  design <- stats::model.matrix(ols)
  ols_beta <- stats::coef(ols)
  ols_sigma2 <- mean(stats::resid(ols)^2)
  starts <- make_start_candidates(ols_beta, ols_sigma2, step_scale)
  score_function <- function(par) semi_score_mean(par, design, response)

  candidates <- lapply(seq_along(starts), function(start_id) {
    root <- tryCatch({
      invisible(utils::capture.output(
        root_result <- suppressWarnings(rootSolve::multiroot(
          f = score_function,
          start = starts[[start_id]],
          maxiter = maxiter,
          rtol = 1e-9,
          atol = 1e-9
        ))
      ))
      root_result$root
    }, error = function(e) NULL)
    if (is.null(root) || any(!is.finite(root))) {
      return(NULL)
    }

    score <- tryCatch(semi_score_i(root, design, response), error = function(e) NULL)
    if (is.null(score) || max(abs(colMeans(score))) > score_tol) {
      return(NULL)
    }

    bread <- tryCatch(
      numeric_bread(root, score_function),
      error = function(e) NULL
    )
    if (is.null(bread)) {
      return(NULL)
    }
    bread_inverse <- tryCatch(
      safe_inverse(bread, rcond_tol = rcond_tol),
      error = function(e) NULL
    )
    if (is.null(bread_inverse)) {
      return(NULL)
    }

    meat <- crossprod(score) / nrow(score)
    covariance <- bread_inverse %*% meat %*% t(bread_inverse) / nrow(score)
    if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
      return(NULL)
    }

    list(
      start_id = start_id,
      par = root,
      score = score,
      bread = bread,
      covariance = covariance,
      distance = root_distance_from_ols(root, ols_beta, ols_sigma2, design)
    )
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) {
    stop("No numerically valid semiparametric root was found.")
  }

  root_key <- vapply(
    candidates,
    function(x) paste(round(x$par, 7), collapse = ":"),
    character(1)
  )
  candidates <- candidates[!duplicated(root_key)]
  selected <- candidates[[which.min(vapply(candidates, `[[`, numeric(1), "distance"))]]
  p <- ncol(design)
  parameter_names <- c(colnames(design), "log_sigma2")
  names(selected$par) <- parameter_names
  rownames(selected$bread) <- colnames(selected$bread) <- parameter_names
  rownames(selected$covariance) <- colnames(selected$covariance) <- parameter_names

  list(
    method = "Semiparametric",
    formula = formula,
    coefficients = stats::setNames(selected$par[seq_len(p)], colnames(design)),
    sigma2 = exp(selected$par[p + 1L]),
    par = selected$par,
    score = selected$score,
    bread = selected$bread,
    covariance = selected$covariance,
    selected_start = selected$start_id,
    candidate_count = length(candidates),
    score_norm = max(abs(colMeans(selected$score))),
    root_distance = selected$distance,
    design = design,
    response = response
  )
}

fit_ols_lm <- function(formula, data, rcond_tol = 1e-10) {
  fit <- stats::lm(formula, data = data)
  design <- stats::model.matrix(fit)
  response <- stats::model.response(stats::model.frame(fit))
  residual <- stats::resid(fit)
  score <- design * residual
  bread <- -crossprod(design) / nrow(design)
  bread_inverse <- safe_inverse(bread, rcond_tol = rcond_tol)
  meat <- crossprod(score) / nrow(design)
  covariance <- bread_inverse %*% meat %*% t(bread_inverse) / nrow(design)

  list(
    method = "OLS",
    formula = formula,
    fit = fit,
    coefficients = stats::coef(fit),
    par = stats::coef(fit),
    score = score,
    bread = bread,
    covariance = covariance,
    design = design,
    response = response
  )
}

stack_mediation_fits <- function(mediator_fit, outcome_fit, baseline = NULL) {
  if (nrow(mediator_fit$score) != nrow(outcome_fit$score)) {
    stop("Mediator and outcome fits must use the same observations.")
  }
  n <- nrow(mediator_fit$score)
  med_names <- paste0("m:", names(mediator_fit$par))
  out_names <- paste0("y:", names(outcome_fit$par))
  parameter <- c(mediator_fit$par, outcome_fit$par)
  names(parameter) <- c(med_names, out_names)
  score <- cbind(mediator_fit$score, outcome_fit$score)
  colnames(score) <- names(parameter)
  bread <- as.matrix(Matrix::bdiag(mediator_fit$bread, outcome_fit$bread))

  if (!is.null(baseline)) {
    baseline <- as.data.frame(baseline)
    if (nrow(baseline) != n) {
      stop("Baseline data must have the same number of rows as the fitted models.")
    }
    baseline_mean <- colMeans(baseline)
    baseline_score <- sweep(as.matrix(baseline), 2, baseline_mean, "-")
    baseline_names <- paste0("xbar:", colnames(baseline))
    parameter <- c(parameter, stats::setNames(baseline_mean, baseline_names))
    score <- cbind(score, baseline_score)
    colnames(score) <- names(parameter)
    bread <- as.matrix(Matrix::bdiag(bread, -diag(ncol(baseline))))
  }

  rownames(bread) <- colnames(bread) <- names(parameter)
  bread_inverse <- safe_inverse(bread)
  meat <- crossprod(score) / n
  covariance <- bread_inverse %*% meat %*% t(bread_inverse) / n
  rownames(covariance) <- colnames(covariance) <- names(parameter)

  list(
    method = mediator_fit$method,
    parameter = parameter,
    score = score,
    bread = bread,
    covariance = covariance,
    mediator_fit = mediator_fit,
    outcome_fit = outcome_fit
  )
}

mediation_effects <- function(stacked_fit, baseline_names = character()) {
  theta <- stacked_fit$parameter
  med_coef_names <- names(stacked_fit$mediator_fit$coefficients)
  out_coef_names <- names(stacked_fit$outcome_fit$coefficients)
  idx <- function(prefix, name) match(paste0(prefix, name), names(theta))

  alpha2 <- unname(theta[idx("m:", "(Intercept)")])
  beta2 <- unname(theta[idx("m:", "T")])
  beta3 <- unname(theta[idx("y:", "T")])
  gamma <- unname(theta[idx("y:", "M")])
  eta <- unname(theta[idx("y:", "T:M")])

  mu0 <- alpha2
  if (length(baseline_names)) {
    med_baseline <- vapply(baseline_names, function(x) {
      unname(theta[idx("m:", x)])
    }, numeric(1))
    baseline_mean <- unname(theta[paste0("xbar:", baseline_names)])
    mu0 <- mu0 + sum(med_baseline * baseline_mean)
  }
  mu1 <- mu0 + beta2

  estimate <- c(
    PNIE = beta2 * gamma,
    TNIE = beta2 * (gamma + eta),
    PNDE = beta3 + eta * mu0,
    TNDE = beta3 + eta * mu1,
    TE = beta3 + beta2 * gamma + eta * mu1
  )

  gradient <- matrix(0, nrow = length(estimate), ncol = length(theta))
  rownames(gradient) <- names(estimate)
  colnames(gradient) <- names(theta)

  gradient["PNIE", idx("m:", "T")] <- gamma
  gradient["PNIE", idx("y:", "M")] <- beta2

  gradient["TNIE", idx("m:", "T")] <- gamma + eta
  gradient["TNIE", idx("y:", "M")] <- beta2
  gradient["TNIE", idx("y:", "T:M")] <- beta2

  gradient["PNDE", idx("m:", "(Intercept)")] <- eta
  gradient["PNDE", idx("y:", "T")] <- 1
  gradient["PNDE", idx("y:", "T:M")] <- mu0

  gradient["TNDE", idx("m:", "(Intercept)")] <- eta
  gradient["TNDE", idx("m:", "T")] <- eta
  gradient["TNDE", idx("y:", "T")] <- 1
  gradient["TNDE", idx("y:", "T:M")] <- mu1

  gradient["TE", idx("m:", "(Intercept)")] <- eta
  gradient["TE", idx("m:", "T")] <- gamma + eta
  gradient["TE", idx("y:", "T")] <- 1
  gradient["TE", idx("y:", "M")] <- beta2
  gradient["TE", idx("y:", "T:M")] <- mu1

  if (length(baseline_names)) {
    for (name in baseline_names) {
      med_index <- idx("m:", name)
      mean_index <- match(paste0("xbar:", name), names(theta))
      med_coefficient <- unname(theta[med_index])
      covariate_mean <- unname(theta[mean_index])
      for (effect in c("PNDE", "TNDE", "TE")) {
        gradient[effect, med_index] <- eta * covariate_mean
        gradient[effect, mean_index] <- eta * med_coefficient
      }
    }
  }

  variance <- gradient %*% stacked_fit$covariance %*% t(gradient)
  std_error <- sqrt(pmax(diag(variance), 0))
  data.frame(
    Effect = names(estimate),
    Estimate = unname(estimate),
    StdError = unname(std_error),
    Lower = unname(estimate - stats::qnorm(0.975) * std_error),
    Upper = unname(estimate + stats::qnorm(0.975) * std_error),
    Method = stacked_fit$method,
    row.names = NULL
  )
}

fit_stacked_mediation <- function(
  data,
  method = c("Semiparametric", "OLS"),
  baseline_names = character(),
  step_scale = 0.10,
  maxiter = 200
) {
  method <- match.arg(method)
  baseline_term <- if (length(baseline_names)) {
    paste(baseline_names, collapse = " + ")
  } else {
    "1"
  }
  mediator_formula <- stats::as.formula(paste("M ~ T +", baseline_term))
  outcome_formula <- stats::as.formula(paste("Y ~ T * M +", baseline_term))
  fit_function <- if (method == "OLS") fit_ols_lm else fit_semiparametric_lm
  fit_arguments <- if (method == "OLS") {
    list()
  } else {
    list(step_scale = step_scale, maxiter = maxiter)
  }

  mediator_fit <- do.call(fit_function, c(list(mediator_formula, data), fit_arguments))
  outcome_fit <- do.call(fit_function, c(list(outcome_formula, data), fit_arguments))
  baseline <- if (length(baseline_names)) data[baseline_names] else NULL
  stacked <- stack_mediation_fits(mediator_fit, outcome_fit, baseline)
  mediation_effects(stacked, baseline_names = baseline_names)
}

generate_standardized_error <- function(
  n,
  scenario = c("gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture")
) {
  scenario <- match.arg(scenario)
  if (scenario == "gaussian") {
    return(stats::rnorm(n))
  }
  if (scenario == "skew_normal") {
    alpha <- 6
    delta <- alpha / sqrt(1 + alpha^2)
    mean_raw <- delta * sqrt(2 / pi)
    sd_raw <- sqrt(1 - 2 * delta^2 / pi)
    raw <- sn::rsn(n, xi = 0, omega = 1, alpha = alpha)
    return((raw - mean_raw) / sd_raw)
  }

  if (scenario == "asymmetric_mixture") {
    probability <- c(0.82, 0.18)
    means <- c(-0.7, 3.2)
    sds <- c(0.35, 0.75)
  } else {
    probability <- c(0.5, 0.5)
    means <- c(-1.7, 1.7)
    sds <- c(0.45, 0.45)
  }
  component <- sample.int(2L, n, replace = TRUE, prob = probability)
  raw <- stats::rnorm(n, mean = means[component], sd = sds[component])
  mixture_mean <- sum(probability * means)
  mixture_variance <- sum(
    probability * (sds^2 + (means - mixture_mean)^2)
  )
  (raw - mixture_mean) / sqrt(mixture_variance)
}

generate_mediation_data <- function(
  n,
  scenario = "gaussian",
  alpha2 = 0.2,
  beta2 = 0.4,
  alpha3 = 0,
  beta3 = 0.5,
  gamma = -0.8,
  eta = 1,
  xi2 = 0.3,
  xi3 = 0.4
) {
  required_cs_packages()
  if (!requireNamespace("sn", quietly = TRUE)) {
    stop("The sn package is required for the skew-normal scenario.")
  }
  X <- stats::rnorm(n)
  T <- stats::rbinom(n, 1, 0.5)
  epsilon_m <- generate_standardized_error(n, scenario)
  epsilon_y <- generate_standardized_error(n, scenario)
  M <- alpha2 + beta2 * T + xi2 * X + epsilon_m
  Y <- alpha3 + beta3 * T + gamma * M + eta * T * M + xi3 * X + epsilon_y
  data.frame(Y = Y, M = M, T = T, X = X)
}

true_mediation_effects <- function(
  alpha2 = 0.2,
  beta2 = 0.4,
  beta3 = 0.5,
  gamma = -0.8,
  eta = 1,
  xi2 = 0.3,
  mean_x = 0
) {
  mu0 <- alpha2 + xi2 * mean_x
  mu1 <- mu0 + beta2
  c(
    PNIE = beta2 * gamma,
    TNIE = beta2 * (gamma + eta),
    PNDE = beta3 + eta * mu0,
    TNDE = beta3 + eta * mu1,
    TE = beta3 + beta2 * gamma + eta * mu1
  )
}

fit_mediation_package <- function(
  data,
  boot = FALSE,
  sims = 499,
  seed = NULL,
  baseline_names = "X"
) {
  required_cs_packages(include_comparators = TRUE)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  baseline_term <- paste(baseline_names, collapse = " + ")
  mediator_formula <- stats::as.formula(paste("M ~ T +", baseline_term))
  outcome_formula <- stats::as.formula(paste("Y ~ T * M +", baseline_term))
  mediator_model <- stats::lm(
    mediator_formula,
    data = data
  )
  outcome_model <- stats::lm(
    outcome_formula,
    data = data
  )
  mediator_model$call$formula <- mediator_formula
  outcome_model$call$formula <- outcome_formula
  fit <- mediation::mediate(
    mediator_model,
    outcome_model,
    treat = "T",
    mediator = "M",
    sims = sims,
    boot = boot,
    boot.ci.type = "perc",
    long = FALSE
  )

  estimate <- c(PNIE = fit$d0, TNIE = fit$d1, PNDE = fit$z0, TNDE = fit$z1, TE = fit$tau.coef)
  intervals <- rbind(fit$d0.ci, fit$d1.ci, fit$z0.ci, fit$z1.ci, fit$tau.ci)
  method <- if (boot) "Imai bootstrap" else "Imai quasi-Bayes"
  data.frame(
    Effect = names(estimate),
    Estimate = unname(estimate),
    StdError = NA_real_,
    Lower = intervals[, 1],
    Upper = intervals[, 2],
    Method = method,
    row.names = NULL
  )
}

fit_medflex <- function(data, baseline_names = "X") {
  required_cs_packages(include_comparators = TRUE)
  centered <- data
  for (name in baseline_names) {
    centered[[name]] <- centered[[name]] - mean(centered[[name]])
  }
  baseline_term <- paste(baseline_names, collapse = " + ")
  imputation_formula <- stats::as.formula(paste("Y ~ T * M +", baseline_term))
  natural_effect_formula <- stats::as.formula(
    paste("Y ~ T0 * T1 +", baseline_term)
  )
  # medflex re-evaluates stored calls in the formula environment.
  evaluation_environment <- new.env(parent = .GlobalEnv)
  evaluation_environment$analysis_data <- centered
  environment(imputation_formula) <- evaluation_environment
  imputation_call <- substitute(
    medflex::neImpute(
      FORMULA,
      family = stats::gaussian(),
      data = analysis_data
    ),
    list(FORMULA = imputation_formula)
  )
  imputed <- eval(imputation_call, envir = evaluation_environment)
  evaluation_environment$expanded_data <- imputed
  evaluation_environment$expData <- imputed
  environment(natural_effect_formula) <- evaluation_environment
  model_call <- substitute(
    medflex::neModel(
      FORMULA,
      family = stats::gaussian(),
      expData = expanded_data,
      se = "robust",
      progress = FALSE
    ),
    list(FORMULA = natural_effect_formula)
  )
  fit <- eval(model_call, envir = evaluation_environment)
  decomposition <- medflex::neEffdecomp(fit)
  test <- summary(decomposition)$test
  source_names <- c(
    PNDE = "pure direct effect",
    TNDE = "total direct effect",
    PNIE = "pure indirect effect",
    TNIE = "total indirect effect",
    TE = "total effect"
  )
  estimate <- test$coefficients[source_names]
  std_error <- test$sigma[source_names]
  estimate <- estimate[c("pure indirect effect", "total indirect effect", "pure direct effect", "total direct effect", "total effect")]
  std_error <- std_error[names(estimate)]
  names(estimate) <- c("PNIE", "TNIE", "PNDE", "TNDE", "TE")
  names(std_error) <- names(estimate)
  data.frame(
    Effect = names(estimate),
    Estimate = unname(estimate),
    StdError = unname(std_error),
    Lower = unname(estimate - stats::qnorm(0.975) * std_error),
    Upper = unname(estimate + stats::qnorm(0.975) * std_error),
    Method = "medflex NEM",
    row.names = NULL
  )
}

scenario_label <- function(scenario) {
  c(
    gaussian = "Gaussian",
    skew_normal = "Skew-normal",
    asymmetric_mixture = "Asymmetric mixture",
    bimodal_mixture = "Symmetric bimodal"
  )[[scenario]]
}

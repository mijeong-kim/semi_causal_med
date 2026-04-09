source("R/jobs_semipara_delta_revision.R")
library(Matrix)
library(sn)

standardize_error_jci <- function(x, mean_x, sd_x) {
  (x - mean_x) / sd_x
}

mixture_moments_jci <- function(prob, means, sds) {
  mix_mean <- sum(prob * means)
  mix_var <- sum(prob * (sds^2 + (means - mix_mean)^2))
  list(mean = mix_mean, sd = sqrt(mix_var))
}

scenario_label_jci <- function(scenario) {
  switch(
    scenario,
    gaussian = "Gaussian",
    skewed = "Skew-normal",
    asymmetric = "Asymmetric mixture",
    bimodal = "Symmetric bimodal",
    stop("Unknown scenario: ", scenario)
  )
}

generate_error_jci <- function(
  n,
  scenario = c("gaussian", "skewed", "asymmetric", "bimodal")
) {
  scenario <- match.arg(scenario)

  if (scenario == "gaussian") {
    e <- stats::rnorm(n)
  } else if (scenario == "skewed") {
    alpha <- 6
    delta <- alpha / sqrt(1 + alpha^2)
    mean_e <- delta * sqrt(2 / pi)
    sd_e <- sqrt(1 - 2 * delta^2 / pi)
    e_raw <- sn::rsn(n, xi = 0, omega = 1, alpha = alpha)
    e <- standardize_error_jci(e_raw, mean_e, sd_e)
  } else if (scenario == "asymmetric") {
    prob <- c(0.82, 0.18)
    means <- c(-0.7, 3.2)
    sds <- c(0.35, 0.75)
    comp <- stats::rbinom(n, size = 1, prob = prob[2]) + 1L
    e_raw <- stats::rnorm(n, mean = means[comp], sd = sds[comp])
    mix_mom <- mixture_moments_jci(prob, means, sds)
    e <- standardize_error_jci(e_raw, mix_mom$mean, mix_mom$sd)
  } else {
    comp <- stats::rbinom(n, size = 1, prob = 0.5)
    means <- c(-1.7, 1.7)
    sds <- c(0.45, 0.45)
    e_raw <- stats::rnorm(n, mean = means[comp + 1L], sd = sds[comp + 1L])
    mix_mom <- mixture_moments_jci(c(0.5, 0.5), means, sds)
    e <- standardize_error_jci(e_raw, mix_mom$mean, mix_mom$sd)
  }
  e
}

generate_interaction_mediation_data <- function(
  n,
  scenario = c("gaussian", "skewed", "asymmetric", "bimodal"),
  alpha2 = 0.2,
  beta2 = 0.4,
  alpha3 = 0.0,
  beta3 = 0.5,
  gamma = -0.8,
  eta = 1.0
) {
  scenario <- match.arg(scenario)
  T <- stats::rbinom(n, size = 1, prob = 0.5)
  e2 <- generate_error_jci(n, scenario = scenario)
  e3 <- generate_error_jci(n, scenario = scenario)

  M <- alpha2 + beta2 * T + e2
  Y <- alpha3 + beta3 * T + gamma * M + eta * T * M + e3

  data.frame(Y = Y, M = M, T = T)
}

interaction_truth <- function(
  alpha2 = 0.2,
  beta2 = 0.4,
  beta3 = 0.5,
  gamma = -0.8,
  eta = 1.0
) {
  c(
    "ACME(0)" = beta2 * gamma,
    "ACME(1)" = beta2 * (gamma + eta),
    "ADE(0)" = beta3 + eta * alpha2,
    "ADE(1)" = beta3 + eta * (alpha2 + beta2),
    "ATE" = beta3 + beta2 * gamma + eta * (alpha2 + beta2)
  )
}

make_zero_padded_jacobian <- function(effect_names, theta) {
  # Build the Jacobian directly on the full stacked parameter vector and
  # leave nuisance coordinates at zero unless an effect depends on them.
  matrix(0, nrow = length(effect_names), ncol = length(theta), dimnames = list(effect_names, names(theta)))
}

interaction_effect_table_simple <- function(
  med_obj,
  out_obj,
  cov_mat,
  include_sigma = FALSE
) {
  med_names <- names(med_obj$coefficients)
  out_names <- names(out_obj$coefficients)
  n_med <- length(med_names)
  shift <- if (include_sigma) n_med + 1L else n_med

  idxm <- function(nm) which(med_names == nm)[1]
  idxo <- function(nm) shift + which(out_names == nm)[1]

  theta <- if (include_sigma) {
    c(med_obj$coefficients, med_obj$sigma2, out_obj$coefficients, out_obj$sigma2)
  } else {
    c(med_obj$coefficients, out_obj$coefficients)
  }

  alpha2 <- unname(theta[idxm("(Intercept)")])
  beta2 <- unname(theta[idxm("T")])
  beta3 <- unname(theta[idxo("T")])
  gamma <- unname(theta[idxo("M")])
  eta <- unname(theta[idxo("T:M")])

  est <- c(
    "ACME(0)" = beta2 * gamma,
    "ACME(1)" = beta2 * (gamma + eta),
    "ADE(0)" = beta3 + eta * alpha2,
    "ADE(1)" = beta3 + eta * (alpha2 + beta2),
    "ATE" = beta3 + beta2 * gamma + eta * (alpha2 + beta2)
  )

  G <- make_zero_padded_jacobian(names(est), theta)
  G["ACME(0)", idxm("T")] <- gamma
  G["ACME(0)", idxo("M")] <- beta2

  G["ACME(1)", idxm("T")] <- gamma + eta
  G["ACME(1)", idxo("M")] <- beta2
  G["ACME(1)", idxo("T:M")] <- beta2

  G["ADE(0)", idxm("(Intercept)")] <- eta
  G["ADE(0)", idxo("T")] <- 1
  G["ADE(0)", idxo("T:M")] <- alpha2

  G["ADE(1)", idxm("(Intercept)")] <- eta
  G["ADE(1)", idxm("T")] <- eta
  G["ADE(1)", idxo("T")] <- 1
  G["ADE(1)", idxo("T:M")] <- alpha2 + beta2

  G["ATE", idxm("(Intercept)")] <- eta
  G["ATE", idxm("T")] <- gamma + eta
  G["ATE", idxo("T")] <- 1
  G["ATE", idxo("M")] <- beta2
  G["ATE", idxo("T:M")] <- alpha2 + beta2

  se <- sqrt(diag(G %*% cov_mat %*% t(G)))

  data.frame(
    Effect = names(est),
    Estimate = unname(est),
    Std.Error = unname(se),
    Lower = unname(est - 1.96 * se),
    Upper = unname(est + 1.96 * se),
    row.names = NULL
  )
}

interaction_effect_table_jobs <- function(
  med_obj,
  out_obj,
  cov_mat,
  xbar,
  include_sigma = FALSE
) {
  med_names <- names(med_obj$coefficients)
  out_names <- names(out_obj$coefficients)
  n_med <- length(med_names)
  shift <- if (include_sigma) n_med + 1L else n_med

  idxm <- function(nm) which(med_names == nm)[1]
  idxo <- function(nm) shift + which(out_names == nm)[1]

  theta <- if (include_sigma) {
    c(med_obj$coefficients, med_obj$sigma2, out_obj$coefficients, out_obj$sigma2)
  } else {
    c(med_obj$coefficients, out_obj$coefficients)
  }

  alpha2 <- unname(theta[idxm("(Intercept)")])
  beta2 <- unname(theta[idxm("treat")])
  xi_eh <- unname(theta[idxm("econ_hard")])
  xi_sex <- unname(theta[idxm("sex")])
  xi_age <- unname(theta[idxm("age")])
  beta3 <- unname(theta[idxo("treat")])
  gamma <- unname(theta[idxo("job_seek")])
  eta <- unname(theta[idxo("treat:job_seek")])
  mu0 <- alpha2 + xi_eh * xbar[1] + xi_sex * xbar[2] + xi_age * xbar[3]
  mu1 <- mu0 + beta2

  est <- c(
    "ACME(0)" = beta2 * gamma,
    "ACME(1)" = beta2 * (gamma + eta),
    "ADE(0)" = beta3 + eta * mu0,
    "ADE(1)" = beta3 + eta * mu1,
    "ATE" = beta3 + beta2 * gamma + eta * mu1
  )

  G <- make_zero_padded_jacobian(names(est), theta)
  G["ACME(0)", idxm("treat")] <- gamma
  G["ACME(0)", idxo("job_seek")] <- beta2

  G["ACME(1)", idxm("treat")] <- gamma + eta
  G["ACME(1)", idxo("job_seek")] <- beta2
  G["ACME(1)", idxo("treat:job_seek")] <- beta2

  G["ADE(0)", idxm("(Intercept)")] <- eta
  G["ADE(0)", idxm("econ_hard")] <- eta * xbar[1]
  G["ADE(0)", idxm("sex")] <- eta * xbar[2]
  G["ADE(0)", idxm("age")] <- eta * xbar[3]
  G["ADE(0)", idxo("treat")] <- 1
  G["ADE(0)", idxo("treat:job_seek")] <- mu0

  G["ADE(1)", idxm("(Intercept)")] <- eta
  G["ADE(1)", idxm("treat")] <- eta
  G["ADE(1)", idxm("econ_hard")] <- eta * xbar[1]
  G["ADE(1)", idxm("sex")] <- eta * xbar[2]
  G["ADE(1)", idxm("age")] <- eta * xbar[3]
  G["ADE(1)", idxo("treat")] <- 1
  G["ADE(1)", idxo("treat:job_seek")] <- mu1

  G["ATE", idxm("(Intercept)")] <- eta
  G["ATE", idxm("treat")] <- gamma + eta
  G["ATE", idxm("econ_hard")] <- eta * xbar[1]
  G["ATE", idxm("sex")] <- eta * xbar[2]
  G["ATE", idxm("age")] <- eta * xbar[3]
  G["ATE", idxo("treat")] <- 1
  G["ATE", idxo("job_seek")] <- beta2
  G["ATE", idxo("treat:job_seek")] <- mu1

  se <- sqrt(diag(G %*% cov_mat %*% t(G)))

  data.frame(
    Effect = names(est),
    Estimate = unname(est),
    Std.Error = unname(se),
    Lower = unname(est - 1.96 * se),
    Upper = unname(est + 1.96 * se),
    row.names = NULL
  )
}

run_jobs_interaction_revision <- function(maxiter = 200, inv_method = c("strict", "ginv")) {
  inv_method <- match.arg(inv_method)
  data(jobs, package = "mediation")
  jobs_data <- jobs[, c("treat", "econ_hard", "sex", "age", "job_seek", "depress2")]
  n <- nrow(jobs_data)
  xbar <- c(
    mean(jobs_data$econ_hard),
    mean(jobs_data$sex),
    mean(jobs_data$age)
  )

  ols_med <- fit_ols_lm_with_score(job_seek ~ treat + econ_hard + sex + age, jobs_data, inv_method = inv_method)
  ols_out <- fit_ols_lm_with_score(depress2 ~ treat * job_seek + econ_hard + sex + age, jobs_data, inv_method = inv_method)
  ols_score <- rbind(ols_med$score, ols_out$score)
  ols_bread <- as.matrix(Matrix::bdiag(ols_med$bread, ols_out$bread))
  ols_bread_inv <- safe_solve(ols_bread, "Stacked OLS Jacobian", inv_method = inv_method)
  ols_cov <- ols_bread_inv %*% (ols_score %*% t(ols_score) / n) %*% t(ols_bread_inv) / n

  semi_med <- fit_semi_indep_lm_stable(job_seek ~ treat + econ_hard + sex + age, jobs_data, maxiter = maxiter, inv_method = inv_method)
  semi_out <- fit_semi_indep_lm_stable(depress2 ~ treat * job_seek + econ_hard + sex + age, jobs_data, maxiter = maxiter, inv_method = inv_method)
  semi_score <- rbind(semi_med$score, semi_out$score)
  semi_bread <- as.matrix(Matrix::bdiag(semi_med$bread, semi_out$bread))
  semi_bread_inv <- safe_solve(semi_bread, "Stacked semiparametric Jacobian", inv_method = inv_method)
  semi_cov <- semi_bread_inv %*% (semi_score %*% t(semi_score) / n) %*% t(semi_bread_inv) / n

  list(
    ols = list(
      med_table = as.data.frame(summary(ols_med$fit)$coef),
      out_table = as.data.frame(summary(ols_out$fit)$coef),
      effects = interaction_effect_table_jobs(ols_med, ols_out, ols_cov, xbar = xbar, include_sigma = FALSE)
    ),
    semi_indep = list(
      med_table = data.frame(
        Estimate = semi_med$coefficients,
        Std.Error = sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)],
        z.value = semi_med$coefficients / sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)],
        p.value = 2 * (1 - pnorm(abs(semi_med$coefficients / sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)])))
      ),
      out_table = data.frame(
        Estimate = semi_out$coefficients,
        Std.Error = sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)],
        z.value = semi_out$coefficients / sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)],
        p.value = 2 * (1 - pnorm(abs(semi_out$coefficients / sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)])))
      ),
      effects = interaction_effect_table_jobs(semi_med, semi_out, semi_cov, xbar = xbar, include_sigma = TRUE)
    )
  )
}

run_uis_interaction_revision <- function(maxiter = 200, inv_method = c("strict", "ginv"), scale_time = 100) {
  inv_method <- match.arg(inv_method)
  data(uis, package = "quantreg")
  raw_data <- stats::na.omit(uis[, c("TREAT", "FRAC", "TIME")])
  uis_data <- data.frame(
    T = raw_data$TREAT,
    M = raw_data$FRAC,
    Y = raw_data$TIME / scale_time
  )
  n <- nrow(uis_data)

  ols_med <- fit_ols_lm_with_score(M ~ T, uis_data, inv_method = inv_method)
  ols_out <- fit_ols_lm_with_score(Y ~ T * M, uis_data, inv_method = inv_method)
  ols_score <- rbind(ols_med$score, ols_out$score)
  ols_bread <- as.matrix(Matrix::bdiag(ols_med$bread, ols_out$bread))
  ols_bread_inv <- safe_solve(ols_bread, "Stacked OLS Jacobian", inv_method = inv_method)
  ols_cov <- ols_bread_inv %*% (ols_score %*% t(ols_score) / n) %*% t(ols_bread_inv) / n

  semi_med <- fit_semi_indep_lm_stable(M ~ T, uis_data, maxiter = maxiter, inv_method = inv_method)
  semi_out <- fit_semi_indep_lm_stable(Y ~ T * M, uis_data, maxiter = maxiter, inv_method = inv_method)
  semi_score <- rbind(semi_med$score, semi_out$score)
  semi_bread <- as.matrix(Matrix::bdiag(semi_med$bread, semi_out$bread))
  semi_bread_inv <- safe_solve(semi_bread, "Stacked semiparametric Jacobian", inv_method = inv_method)
  semi_cov <- semi_bread_inv %*% (semi_score %*% t(semi_score) / n) %*% t(semi_bread_inv) / n

  list(
    scale_time = scale_time,
    ols = list(
      med_table = as.data.frame(summary(ols_med$fit)$coef),
      out_table = as.data.frame(summary(ols_out$fit)$coef),
      effects = interaction_effect_table_simple(ols_med, ols_out, ols_cov, include_sigma = FALSE)
    ),
    semi_indep = list(
      med_table = data.frame(
        Estimate = semi_med$coefficients,
        Std.Error = sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)],
        z.value = semi_med$coefficients / sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)],
        p.value = 2 * (1 - pnorm(abs(semi_med$coefficients / sqrt(diag(semi_med$covariance))[seq_along(semi_med$coefficients)])))
      ),
      out_table = data.frame(
        Estimate = semi_out$coefficients,
        Std.Error = sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)],
        z.value = semi_out$coefficients / sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)],
        p.value = 2 * (1 - pnorm(abs(semi_out$coefficients / sqrt(diag(semi_out$covariance))[seq_along(semi_out$coefficients)])))
      ),
      effects = interaction_effect_table_simple(semi_med, semi_out, semi_cov, include_sigma = TRUE)
    )
  )
}

summarize_simulation_jci <- function(records, true_effects) {
  within_tol <- function(x, target) abs(x - target)

  do.call(
    rbind,
    lapply(split(records, list(records$Scenario, records$Method, records$Effect), drop = TRUE), function(df) {
      truth <- true_effects[df$Effect[1]]
      data.frame(
        Scenario = df$Scenario[1],
        Method = df$Method[1],
        Effect = df$Effect[1],
        Bias = mean(df$Estimate - truth),
        RMSE = sqrt(mean((df$Estimate - truth)^2)),
        Coverage95 = mean(df$Lower <= truth & df$Upper >= truth),
        AvgLength = mean(df$Upper - df$Lower),
        SuccessRate = unique(df$SuccessRate)[1],
        stringsAsFactors = FALSE
      )
    })
  )
}

run_interaction_simulation <- function(
  n = 300,
  nsim = 1000,
  scenarios = c("gaussian", "skewed", "asymmetric", "bimodal"),
  alpha2 = 0.2,
  beta2 = 0.4,
  alpha3 = 0.0,
  beta3 = 0.5,
  gamma = -0.8,
  eta = 1.0,
  maxiter = 150,
  seed = 20260329
) {
  set.seed(seed)
  truth <- interaction_truth(alpha2 = alpha2, beta2 = beta2, beta3 = beta3, gamma = gamma, eta = eta)

  records <- list()
  counter <- 1L

  for (scenario in scenarios) {
    success_counts <- c(OLS = 0L, Semiparametric = 0L)

    for (s in seq_len(nsim)) {
      dat <- generate_interaction_mediation_data(
        n = n,
        scenario = scenario,
        alpha2 = alpha2,
        beta2 = beta2,
        alpha3 = alpha3,
        beta3 = beta3,
        gamma = gamma,
        eta = eta
      )

      ols_fit <- tryCatch({
        med <- fit_ols_lm_with_score(M ~ T, dat)
        out <- fit_ols_lm_with_score(Y ~ T * M, dat)
        bread <- as.matrix(Matrix::bdiag(med$bread, out$bread))
        inv_bread <- safe_solve(bread, "Stacked OLS Jacobian")
        score <- rbind(med$score, out$score)
        cov_mat <- inv_bread %*% (score %*% t(score) / n) %*% t(inv_bread) / n
        interaction_effect_table_simple(med, out, cov_mat, include_sigma = FALSE)
      }, error = function(e) NULL)

      if (!is.null(ols_fit)) {
        success_counts["OLS"] <- success_counts["OLS"] + 1L
        for (j in seq_len(nrow(ols_fit))) {
          records[[counter]] <- data.frame(
            Scenario = scenario_label_jci(scenario),
            Method = "OLS",
            Effect = ols_fit$Effect[j],
            Estimate = ols_fit$Estimate[j],
            Lower = ols_fit$Lower[j],
            Upper = ols_fit$Upper[j],
            SuccessRate = NA_real_,
            stringsAsFactors = FALSE
          )
          counter <- counter + 1L
        }
      }

      semi_fit <- tryCatch({
        med <- fit_semi_indep_lm_stable(M ~ T, dat, maxiter = maxiter)
        out <- fit_semi_indep_lm_stable(Y ~ T * M, dat, maxiter = maxiter)
        bread <- as.matrix(Matrix::bdiag(med$bread, out$bread))
        inv_bread <- safe_solve(bread, "Stacked semiparametric Jacobian")
        score <- rbind(med$score, out$score)
        cov_mat <- inv_bread %*% (score %*% t(score) / n) %*% t(inv_bread) / n
        interaction_effect_table_simple(med, out, cov_mat, include_sigma = TRUE)
      }, error = function(e) NULL)

      if (!is.null(semi_fit)) {
        success_counts["Semiparametric"] <- success_counts["Semiparametric"] + 1L
        for (j in seq_len(nrow(semi_fit))) {
          records[[counter]] <- data.frame(
            Scenario = scenario_label_jci(scenario),
            Method = "Semiparametric",
            Effect = semi_fit$Effect[j],
            Estimate = semi_fit$Estimate[j],
            Lower = semi_fit$Lower[j],
            Upper = semi_fit$Upper[j],
            SuccessRate = NA_real_,
            stringsAsFactors = FALSE
          )
          counter <- counter + 1L
        }
      }
    }

    scenario_name <- scenario_label_jci(scenario)
    for (meth in names(success_counts)) {
      idx <- vapply(records, function(x) identical(x$Scenario[1], scenario_name) && identical(x$Method[1], meth), logical(1))
      for (k in which(idx)) {
        records[[k]]$SuccessRate <- success_counts[meth] / nsim
      }
    }
  }

  records_df <- do.call(rbind, records)
  summary_df <- summarize_simulation_jci(records_df, truth)
  rownames(summary_df) <- NULL

  list(
    truth = truth,
    records = records_df,
    summary = summary_df
  )
}

run_boundary_power_example_jci <- function(
  n = 220,
  nsim = 1000,
  scenario = "asymmetric",
  alpha2 = 0.2,
  beta2 = 0.26,
  alpha3 = 0.0,
  beta3 = 0.5,
  gamma = -0.26,
  eta = 0.8,
  effect = "ACME(0)",
  maxiter = 120,
  seed = 20260415
) {
  sim_res <- run_interaction_simulation(
    n = n,
    nsim = nsim,
    scenarios = scenario,
    alpha2 = alpha2,
    beta2 = beta2,
    alpha3 = alpha3,
    beta3 = beta3,
    gamma = gamma,
    eta = eta,
    maxiter = maxiter,
    seed = seed
  )

  effect_records <- subset(sim_res$records, Effect == effect)
  out <- aggregate(
    cbind(
      MeanEstimate = effect_records$Estimate,
      RejectZero = as.numeric(effect_records$Lower > 0 | effect_records$Upper < 0),
      AvgLength = effect_records$Upper - effect_records$Lower
    ) ~ Method,
    data = effect_records,
    FUN = mean
  )

  out$Scenario <- scenario_label_jci(scenario)
  out$Effect <- effect
  out$Truth <- unname(sim_res$truth[effect])
  out$SampleSize <- n
  out$beta2 <- beta2
  out$gamma <- gamma
  out$eta <- eta
  out <- out[, c("Scenario", "SampleSize", "Effect", "Truth", "beta2", "gamma", "eta", "Method", "MeanEstimate", "RejectZero", "AvgLength")]
  rownames(out) <- NULL
  out
}

if (sys.nframe() == 0) {
  jobs_res <- run_jobs_interaction_revision()
  uis_res <- run_uis_interaction_revision()
  sim_res <- run_interaction_simulation(n = 300, nsim = 1000)
  boundary_res <- run_boundary_power_example_jci()

  utils::write.csv(
    jobs_res$ols$effects,
    file = "results/jci_jobs_interaction_ols_effects.csv",
    row.names = FALSE
  )
  utils::write.csv(
    jobs_res$semi_indep$effects,
    file = "results/jci_jobs_interaction_semiparametric_effects.csv",
    row.names = FALSE
  )
  utils::write.csv(
    sim_res$summary,
    file = "results/jci_interaction_simulation_summary.csv",
    row.names = FALSE
  )
  utils::write.csv(
    uis_res$ols$effects,
    file = "results/jci_uis_interaction_ols_effects.csv",
    row.names = FALSE
  )
  utils::write.csv(
    uis_res$semi_indep$effects,
    file = "results/jci_uis_interaction_semiparametric_effects.csv",
    row.names = FALSE
  )
  utils::write.csv(
    boundary_res,
    file = "results/jci_boundary_power_example.csv",
    row.names = FALSE
  )

  cat("\nOLS interaction effects\n")
  print(jobs_res$ols$effects)
  cat("\nSemiparametric interaction effects\n")
  print(jobs_res$semi_indep$effects)
  cat("\nUIS OLS interaction effects\n")
  print(uis_res$ols$effects)
  cat("\nUIS semiparametric interaction effects\n")
  print(uis_res$semi_indep$effects)
  cat("\nInteraction simulation summary\n")
  print(sim_res$summary)
  cat("\nBoundary power example\n")
  print(boundary_res)
}

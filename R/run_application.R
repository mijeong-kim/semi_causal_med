source(file.path("R", "semiparametric_mediation.R"))

required_jkss_packages(include_comparators = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

data_environment <- new.env(parent = emptyenv())
load(file.path("data", "jobs.RData"), envir = data_environment)
jobs <- data_environment$jobs
analysis_data <- stats::na.omit(data.frame(
  Y = jobs$depress2,
  M = jobs$job_seek,
  T = jobs$treat,
  depress1 = jobs$depress1,
  econ_hard = jobs$econ_hard,
  sex = jobs$sex,
  age = jobs$age
))
baseline_names <- c("depress1", "econ_hard", "sex", "age")

fit_method <- function(method, step_scale = 0.10,
                       covariates = baseline_names) {
  mediator_formula <- stats::reformulate(c("T", covariates), response = "M")
  outcome_formula <- stats::reformulate(c("T * M", covariates), response = "Y")
  if (method == "OLS") {
    mediator <- fit_ols_lm(mediator_formula, analysis_data)
    outcome <- fit_ols_lm(outcome_formula, analysis_data)
  } else {
    mediator <- fit_semiparametric_lm(
      mediator_formula, analysis_data, step_scale = step_scale
    )
    outcome <- fit_semiparametric_lm(
      outcome_formula, analysis_data, step_scale = step_scale
    )
  }
  stacked <- stack_mediation_fits(
    mediator,
    outcome,
    analysis_data[covariates]
  )
  list(
    mediator = mediator,
    outcome = outcome,
    effects = mediation_effects(stacked, covariates)
  )
}

message("Fitting the baseline-adjusted JOBS II illustration.")
ols <- fit_method("OLS")
semiparametric <- fit_method("Semiparametric")
set.seed(20260831)
imai_quasi <- fit_mediation_package(
  analysis_data,
  boot = FALSE,
  sims = 5000,
  seed = 20260831,
  baseline_names = baseline_names
)
invisible(utils::capture.output(
  imai_bootstrap <- suppressMessages(fit_mediation_package(
    analysis_data,
    boot = TRUE,
    sims = 1999,
    seed = 20260901,
    baseline_names = baseline_names
  ))
))
medflex <- suppressMessages(fit_medflex(analysis_data, baseline_names))

application_effects <- do.call(rbind, list(
  ols$effects,
  semiparametric$effects,
  imai_quasi,
  imai_bootstrap,
  medflex
))
application_effects$SampleSize <- nrow(analysis_data)
utils::write.csv(
  application_effects,
  "results/application_effects.csv",
  row.names = FALSE
)

coefficient_table <- function(fit, equation) {
  number_of_coefficients <- length(fit$coefficients)
  data.frame(
    Equation = equation,
    Method = fit$method,
    Term = names(fit$coefficients),
    Estimate = unname(fit$coefficients),
    StdError = sqrt(diag(fit$covariance))[seq_len(number_of_coefficients)],
    stringsAsFactors = FALSE
  )
}
application_regression <- do.call(rbind, list(
  coefficient_table(ols$mediator, "Mediator"),
  coefficient_table(semiparametric$mediator, "Mediator"),
  coefficient_table(ols$outcome, "Outcome"),
  coefficient_table(semiparametric$outcome, "Outcome")
))
utils::write.csv(
  application_regression,
  "results/application_regression.csv",
  row.names = FALSE
)

residual_diagnostics <- function(fit, equation) {
  residual <- drop(fit$response - fit$design %*% fit$coefficients)
  standardized <- (residual - mean(residual)) / stats::sd(residual)
  fitted_value <- drop(fit$design %*% fit$coefficients)
  variance_fit <- stats::lm(I(residual^2) ~ fitted_value)
  bp_statistic <- length(residual) * summary(variance_fit)$r.squared
  predictors <- fit$design[, colnames(fit$design) != "(Intercept)", drop = FALSE]
  full_variance_fit <- stats::lm(I(residual^2) ~ predictors)
  full_df <- full_variance_fit$rank - 1L
  full_lm <- length(residual) * summary(full_variance_fit)$r.squared
  data.frame(
    Equation = equation,
    Method = fit$method,
    Mean = mean(residual),
    SD = stats::sd(residual),
    Skewness = mean(standardized^3),
    ExcessKurtosis = mean(standardized^4) - 3,
    ShapiroW = unname(stats::shapiro.test(residual)$statistic),
    ShapiroP = stats::shapiro.test(residual)$p.value,
    BPStatistic = bp_statistic,
    BPP = stats::pchisq(bp_statistic, df = 1, lower.tail = FALSE),
    FullVarianceLM = full_lm,
    FullVarianceDF = full_df,
    FullVarianceP = stats::pchisq(full_lm, df = full_df, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}
application_diagnostics <- rbind(
  residual_diagnostics(ols$mediator, "Mediator"),
  residual_diagnostics(ols$outcome, "Outcome")
)
utils::write.csv(
  application_diagnostics,
  "results/application_diagnostics.csv",
  row.names = FALSE
)

support <- do.call(rbind, lapply(c("M", "Y"), function(variable) {
  x <- analysis_data[[variable]]
  data.frame(Variable = variable, N = length(x), UniqueValues = length(unique(x)),
             Minimum = min(x), Maximum = max(x),
             AtMinimum = mean(x == min(x)), AtMaximum = mean(x == max(x)))
}))
utils::write.csv(support, "results/application_support.csv", row.names = FALSE)

message("Checking the earlier JOBS II specification without baseline depression.")
previous_covariates <- setdiff(baseline_names, "depress1")
previous_ols <- fit_method("OLS", covariates = previous_covariates)
previous_sp <- fit_method("Semiparametric", covariates = previous_covariates)
covariate_sensitivity <- rbind(
  transform(ols$effects, Specification = "With baseline depression"),
  transform(semiparametric$effects, Specification = "With baseline depression"),
  transform(previous_ols$effects, Specification = "Without baseline depression"),
  transform(previous_sp$effects, Specification = "Without baseline depression")
)
utils::write.csv(covariate_sensitivity,
                 "results/application_covariate_sensitivity.csv", row.names = FALSE)

scales <- c(0.05, 0.10, 0.20)
application_sensitivity <- do.call(rbind, lapply(scales, function(scale) {
  fit <- if (scale == 0.10) semiparametric else
    fit_method("Semiparametric", step_scale = scale)
  pnie <- fit$effects[fit$effects$Effect == "PNIE", ]
  data.frame(
    StepScale = scale,
    PNIE = pnie$Estimate,
    Lower = pnie$Lower,
    Upper = pnie$Upper,
    MediatorSelectedStart = fit$mediator$selected_start,
    OutcomeSelectedStart = fit$outcome$selected_start,
    MediatorCandidateRoots = fit$mediator$candidate_count,
    OutcomeCandidateRoots = fit$outcome$candidate_count,
    MediatorScoreNorm = fit$mediator$score_norm,
    OutcomeScoreNorm = fit$outcome$score_norm,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  application_sensitivity,
  "results/application_algorithm_sensitivity.csv",
  row.names = FALSE
)

grDevices::pdf("figures/application_effects.pdf", width = 8.0, height = 5.6)
old_par <- graphics::par(no.readonly = TRUE)
effect_order <- c("TE", "TNDE", "PNDE", "TNIE", "PNIE")
method_order <- c(
  "OLS", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM", "Semiparametric"
)
point_symbols <- c(1, 2, 0, 5, 16)
line_types <- c(1, 2, 3, 4, 1)
# A separate legend panel keeps labels clear of every confidence interval.
graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(0.8, 4.8))
graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::legend(
  "center",
  legend = method_order,
  pch = point_symbols,
  lty = line_types,
  bty = "n",
  ncol = 3,
  cex = 0.92,
  text.width = 1.18 * max(graphics::strwidth(method_order, cex = 0.92)),
  y.intersp = 1.4,
  x.intersp = 1.0
)
graphics::par(mar = c(4.2, 5.2, 0.4, 0.8), las = 1)
plot_data <- application_effects[
  match(
    as.vector(t(outer(effect_order, method_order, paste, sep = "::"))),
    paste(application_effects$Effect, application_effects$Method, sep = "::")
  ),
]
base_y <- rep(seq_along(effect_order), each = length(method_order))
offset <- rep(seq(-0.28, 0.28, length.out = length(method_order)), length(effect_order))
y <- base_y + offset
limits <- range(c(plot_data$Lower, plot_data$Upper), finite = TRUE)
graphics::plot(
  NA,
  xlim = limits,
  ylim = c(0.5, length(effect_order) + 0.5),
  xlab = "Effect on post-intervention depressive symptoms (95% interval)",
  ylab = "",
  yaxt = "n",
  bty = "n"
)
graphics::axis(2, at = seq_along(effect_order), labels = effect_order, tick = FALSE)
graphics::abline(v = 0, lty = 3, col = "grey55")
for (j in seq_along(method_order)) {
  index <- seq(j, nrow(plot_data), by = length(method_order))
  graphics::segments(
    plot_data$Lower[index], y[index], plot_data$Upper[index], y[index],
    lty = line_types[j], col = if (j == 5) "black" else "grey35"
  )
  graphics::points(
    plot_data$Estimate[index], y[index], pch = point_symbols[j],
    col = if (j == 5) "black" else "grey25", cex = 0.9
  )
}
graphics::layout(1)
graphics::par(old_par)
grDevices::dev.off()

grDevices::pdf("figures/application_residuals.pdf", width = 7.6, height = 6.6)
old_par <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 1.6, 0.6))
for (entry in list(
  list(fit = ols$mediator, title = "Mediator model"),
  list(fit = ols$outcome, title = "Outcome model")
)) {
  residual <- stats::resid(entry$fit$fit)
  density <- stats::density(residual)
  grid <- seq(min(density$x), max(density$x), length.out = 300)
  gaussian_density <- stats::dnorm(grid, mean(residual), stats::sd(residual))
  graphics::plot(
    density,
    main = entry$title,
    xlab = "OLS residual",
    ylab = "Density",
    ylim = c(0, 1.30 * max(density$y, gaussian_density)),
    lwd = 1.5
  )
  graphics::lines(
    grid,
    gaussian_density,
    lty = 2,
    lwd = 1.2
  )
  graphics::legend(
    "topright",
    legend = c("Kernel density", "Matched Gaussian"),
    lty = c(1, 2),
    bty = "n",
    cex = 0.8
  )
  fitted_value <- stats::fitted(entry$fit$fit)
  graphics::plot(fitted_value, abs(residual), pch = 16, cex = 0.45,
                 col = grDevices::adjustcolor("grey25", alpha.f = 0.25),
                 xlab = "OLS fitted value", ylab = "Absolute OLS residual",
                 main = paste(entry$title, "spread"))
  graphics::lines(stats::lowess(fitted_value, abs(residual)), lwd = 2)
}
graphics::par(old_par)
grDevices::dev.off()

writeLines(c(
  "Source: jobs data distributed with mediation 4.5.1.",
  paste("Rows loaded:", nrow(jobs), "Complete analysis rows:", nrow(analysis_data)),
  "Treatment: randomized assignment (treat); mediator: job_seek; outcome: depress2.",
  paste("Baseline covariates:", paste(baseline_names, collapse = ", ")),
  "Covariate sensitivity removes baseline depression while retaining the other terms.",
  "Variance diagnostics: n R^2 from squared OLS residuals on fitted values (1 df)",
  "and separately on all non-intercept columns of each original design matrix.",
  "The package illustration contains complete measurements for all selected variables.",
  "All observations, transformations and model specifications were fixed before fitting."
), "results/application_analysis_specification.txt")
writeLines(utils::capture.output(sessionInfo()), "results/application_sessionInfo.txt")

cat("Application results and figures written to results/ and figures/.\n")

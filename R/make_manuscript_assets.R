source(file.path("R", "results_io.R"))
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
source(file.path("R", "semiparametric_mediation.R"))

read_result <- function(name) {
  path <- file.path("results", name)
  if (!jkss_result_exists(path)) {
    stop("Missing result file: ", path)
  }
  read_jkss_csv(path, check.names = FALSE)
}

format_number <- function(x, digits = 3) {
  rounded <- round(x, digits)
  rounded[!is.na(rounded) & abs(rounded) < 0.5 * 10^(-digits)] <- 0
  ifelse(is.na(rounded), "--", formatC(rounded, format = "f", digits = digits))
}

format_bold_number <- function(x, bold, digits = 3) {
  value <- format_number(x, digits)
  ifelse(bold & !is.na(x), paste0("\\textbf{", value, "}"), value)
}

format_interval <- function(estimate, lower, upper) {
  paste0(
    format_number(estimate), " [", format_number(lower), ", ",
    format_number(upper), "]"
  )
}

method_latex <- function(x) {
  labels <- c(
    "OLS" = "OLS delta",
    "Semiparametric" = "Proposed",
    "Imai quasi-Bayes" = "Imai quasi-Bayes",
    "Imai bootstrap" = "Imai bootstrap",
    "medflex NEM" = "\\texttt{medflex} NEM"
  )
  unname(labels[x])
}

method_latex_full <- function(x) {
  labels <- c(
    "OLS" = "OLS",
    "Semiparametric" = "Proposed",
    "Imai quasi-Bayes" = "Imai-QB",
    "Imai bootstrap" = "Imai-Boot",
    "medflex NEM" = "NEM"
  )
  unname(labels[x])
}

scenario_latex_full <- function(x) {
  labels <- c(
    "Gaussian" = "G",
    "Skew-normal" = "SN",
    "Asymmetric mixture" = "AM",
    "Symmetric bimodal" = "BM"
  )
  unname(labels[x])
}

write_lines <- function(lines, path) {
  writeLines(lines, path, useBytes = TRUE)
}

is_minimum <- function(value, candidates) {
  best <- min(candidates, na.rm = TRUE)
  is.finite(value) && abs(value - best) < 1e-12
}

main <- read_result("main_simulation_summary.csv")
comparator <- read_result("comparator_summary.csv")
power <- read_result("size_power_summary.csv")
sensitivity <- read_result("algorithm_sensitivity_summary.csv")
application <- read_result("application_effects.csv")
regression <- read_result("application_regression.csv")
diagnostics <- read_result("application_diagnostics.csv")
application_sensitivity <- read_result("application_algorithm_sensitivity.csv")
application_covariates <- read_result("application_covariate_sensitivity.csv")

make_main_table <- function() {
  selected <- main[main$Effect %in% c("PNIE", "TE"), ]
  keys <- unique(selected[c("SampleSize", "Scenario", "Method")])
  scenario_order <- c(
    "Gaussian", "Skew-normal", "Asymmetric mixture", "Symmetric bimodal"
  )
  method_order <- c("OLS", "Semiparametric")
  keys <- keys[order(
    keys$SampleSize,
    match(keys$Scenario, scenario_order),
    match(keys$Method, method_order)
  ), ]
  rows <- vapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    pnie <- selected[
      selected$SampleSize == key$SampleSize &
        selected$Scenario == key$Scenario &
        selected$Method == key$Method & selected$Effect == "PNIE",
    ]
    te <- selected[
      selected$SampleSize == key$SampleSize &
        selected$Scenario == key$Scenario &
        selected$Method == key$Method & selected$Effect == "TE",
    ]
    pnie_group <- selected[
      selected$SampleSize == key$SampleSize &
        selected$Scenario == key$Scenario & selected$Effect == "PNIE",
    ]
    te_group <- selected[
      selected$SampleSize == key$SampleSize &
        selected$Scenario == key$Scenario & selected$Effect == "TE",
    ]
    paste(
      key$SampleSize,
      scenario_latex_full(key$Scenario),
      method_latex(key$Method),
      format_number(pnie$Bias),
      format_bold_number(
        pnie$RMSE,
        is_minimum(pnie$RMSE, pnie_group$RMSE)
      ),
      format_number(pnie$Coverage95),
      format_bold_number(
        pnie$AvgLength,
        is_minimum(pnie$AvgLength, pnie_group$AvgLength)
      ),
      format_number(te$Bias),
      format_bold_number(
        te$RMSE,
        is_minimum(te$RMSE, te_group$RMSE)
      ),
      format_number(te$Coverage95),
      format_bold_number(
        te$AvgLength,
        is_minimum(te$AvgLength, te_group$AvgLength)
      ),
      sep = " & "
    )
  }, character(1))
  rows <- paste0(rows, " \\\\")
  breaks <- which(head(keys$Scenario, -1) != tail(keys$Scenario, -1) |
                  head(keys$SampleSize, -1) != tail(keys$SampleSize, -1))
  rows[breaks] <- paste0(rows[breaks], "\n\\midrule")
  write_lines(c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Main Monte Carlo results based on 1,000 attempts per cell and summarized over valid fits; counts and success rates appear in Table~S1. G, Gaussian; SN, skew-normal; AM, asymmetric mixture; BM, symmetric bimodal. Boldface marks lower RMSE and shorter length within each design and effect; these descriptive highlights must be interpreted jointly with coverage and numerical success. Bias is signed}",
    "\\label{tab:main_simulation}",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\begin{tabular}{rllrrrrrrrr}",
    "\\toprule",
    "& & & \\multicolumn{4}{c}{PNIE} & \\multicolumn{4}{c}{TE} \\\\",
    "\\cmidrule(lr){4-7}\\cmidrule(lr){8-11}",
    "$n$ & Error & Method & Bias & RMSE & Cov. & Length & Bias & RMSE & Cov. & Length \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/main_simulation_table.tex")
}

make_comparator_table <- function() {
  selected <- comparator[comparator$Effect %in% c("PNIE", "TE"), ]
  keys <- unique(selected[c("Scenario", "Method")])
  scenario_order <- c(
    "Gaussian", "Skew-normal", "Asymmetric mixture", "Symmetric bimodal"
  )
  method_order <- c(
    "OLS", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM", "Semiparametric"
  )
  keys <- keys[order(
    match(keys$Scenario, scenario_order),
    match(keys$Method, method_order)
  ), ]
  rows <- vapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, ]
    pnie <- selected[
      selected$Scenario == key$Scenario & selected$Method == key$Method &
        selected$Effect == "PNIE",
    ]
    te <- selected[
      selected$Scenario == key$Scenario & selected$Method == key$Method &
        selected$Effect == "TE",
    ]
    pnie_group <- selected[
      selected$Scenario == key$Scenario & selected$Effect == "PNIE",
    ]
    te_group <- selected[
      selected$Scenario == key$Scenario & selected$Effect == "TE",
    ]
    paste(
      scenario_latex_full(key$Scenario),
      method_latex(key$Method),
      format_bold_number(
        pnie$RMSE,
        is_minimum(round(pnie$RMSE, 3), round(pnie_group$RMSE, 3))
      ),
      format_number(pnie$Coverage95),
      format_bold_number(
        pnie$AvgLength,
        is_minimum(round(pnie$AvgLength, 3), round(pnie_group$AvgLength, 3))
      ),
      format_bold_number(
        te$RMSE,
        is_minimum(round(te$RMSE, 3), round(te_group$RMSE, 3))
      ),
      format_number(te$Coverage95),
      format_bold_number(
        te$AvgLength,
        is_minimum(round(te$AvgLength, 3), round(te_group$AvgLength, 3))
      ),
      sep = " & "
    )
  }, character(1))
  rows <- paste0(rows, " \\\\")
  breaks <- which(head(keys$Scenario, -1) != tail(keys$Scenario, -1))
  rows[breaks] <- paste0(rows[breaks], "\n\\midrule")
  write_lines(c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Comparator study at $n=300$ based on 1,000 attempts per cell, summarized over valid fits. The Imai procedures use 1,000 quasi-Bayesian draws or 499 bootstrap resamples per data set. Error abbreviations are as in Table~\\ref{tab:main_simulation}. Boldface marks the smallest displayed RMSE and shortest length within each design and effect; these descriptive highlights are interpreted jointly with coverage and numerical success}",
    "\\label{tab:comparators}",
    "\\small",
    "\\begin{tabular}{llrrrrrr}",
    "\\toprule",
    "& & \\multicolumn{3}{c}{PNIE} & \\multicolumn{3}{c}{TE} \\\\",
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
    "Error & Method & RMSE & Cov. & Length & RMSE & Cov. & Length \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/comparator_table.tex")
}

make_application_table <- function() {
  effect_order <- c("PNIE", "TNIE", "PNDE", "TNDE", "TE")
  rows <- vapply(effect_order, function(effect) {
    result <- application[
      application$Method == "Semiparametric" & application$Effect == effect,
    ]
    stopifnot(nrow(result) == 1L)
    paste(
      effect,
      format_number(result$Estimate, 4),
      format_number(result$Lower, 4),
      format_number(result$Upper, 4),
      sep = " & "
    )
  }, character(1))
  rows <- paste0(rows, " \\\\")
  write_lines(c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Proposed density-score estimates of JOBS II natural effects and nominal 95\\% limits ($n=899$). Negative values indicate lower post-intervention depressive symptoms under assignment to the intervention}",
    "\\label{tab:application_effects}",
    "\\small",
    "\\begin{tabular}{lrrr}",
    "\\toprule",
    "Effect & Estimate & Lower & Upper \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/application_effects_table.tex")
}

make_application_full_table <- function() {
  method_order <- c(
    "OLS", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM", "Semiparametric"
  )
  effect_order <- c("PNIE", "TNIE", "PNDE", "TNDE", "TE")
  rows <- unlist(lapply(seq_along(effect_order), function(effect_index) {
    effect <- effect_order[effect_index]
    effect_rows <- vapply(method_order, function(method) {
      subset <- application[
        application$Method == method & application$Effect == effect,
      ]
      paste(
        effect,
        method_latex(method),
        format_number(subset$Estimate, 4),
        format_number(subset$Lower, 4),
        format_number(subset$Upper, 4),
        sep = " & "
      )
    }, character(1))
    effect_rows <- paste0(effect_rows, " \\\\")
    if (effect_index < length(effect_order)) {
      c(effect_rows, "\\addlinespace[2pt]")
    } else {
      effect_rows
    }
  }))
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Complete JOBS II natural-effect estimates and nominal 95\\% limits for five methods ($n=899$). Negative values indicate lower post-intervention depressive symptoms under assignment to the intervention}",
    "\\label{tab:application_effects_full}",
    "\\small",
    "\\begin{tabular}{llrrr}",
    "\\toprule",
    "Effect & Method & Estimate & Lower & Upper \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/application_effects_full_table.tex")
}

make_power_table <- function() {
  ordered <- power[order(power$Beta2, power$Method), ]
  rows <- paste0(vapply(seq_len(nrow(ordered)), function(i) {
    x <- ordered[i, ]
    paste(
      format_number(x$Beta2, 2),
      format_number(x$TruePNIE, 3),
      method_latex(x$Method),
      x$Replications,
      format_number(x$SuccessRate),
      format_number(x$RejectionRate),
      format_number(x$MonteCarloSE),
      sep = " & "
    )
  }, character(1)), " \\\\")
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Empirical size and power for PNIE under the standardized asymmetric-mixture error, $n=300$, based on 1,000 attempted replications. Rejection and MCSE are conditional on numerically successful fits; Success uses all attempts}",
    "\\label{tab:size_power}",
    "\\begin{tabular}{rrlrrrr}",
    "\\toprule",
    "$\\beta_2$ & True PNIE & Method & Reps. & Success & Reject & MCSE \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/size_power_table.tex")
}

make_sensitivity_table <- function() {
  ordered <- sensitivity[order(sensitivity$Scenario, sensitivity$StepScale), ]
  rows <- paste0(vapply(seq_len(nrow(ordered)), function(i) {
    x <- ordered[i, ]
    paste(
      x$Scenario,
      format_number(x$StepScale, 2),
      x$Replications,
      format_number(x$SuccessRate),
      format_number(x$MeanPNIE),
      format_number(x$SDEstimate),
      format_number(x$MeanCandidateRootsMediator, 2),
      format_number(x$MeanCandidateRootsOutcome, 2),
      sep = " & "
    )
  }, character(1)), " \\\\")
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Sensitivity to the coefficient perturbation scale at $n=200$, based on 300 common data sets per error distribution. Success uses all attempts; other summaries use numerically successful fits}",
    "\\label{tab:algorithm_sensitivity}",
    "\\small",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    "Error & Scale & Reps. & Success & Mean PNIE & SD & Roots M & Roots Y \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/algorithm_sensitivity_table.tex")
}

make_regression_table <- function() {
  ordered <- regression[order(
    match(regression$Equation, c("Mediator", "Outcome")),
    match(regression$Method, c("OLS", "Semiparametric"))
  ), ]
  rows <- paste0(vapply(seq_len(nrow(ordered)), function(i) {
    x <- ordered[i, ]
    paste(
      x$Equation,
      method_latex(x$Method),
      gsub("_", "\\_", x$Term, fixed = TRUE),
      format_number(x$Estimate, 4),
      format_number(x$StdError, 4),
      sep = " & "
    )
  }, character(1)), " \\\\")
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Regression estimates for the JOBS II analysis. Variable names follow the reproducibility code}",
    "\\label{tab:application_regression}",
    "\\footnotesize",
    "\\begin{tabular}{lllrr}",
    "\\toprule",
    "Equation & Method & Term & Estimate & SE \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/application_regression_table.tex")
}

make_diagnostics_table <- function() {
  rows <- paste0(vapply(seq_len(nrow(diagnostics)), function(i) {
    x <- diagnostics[i, ]
    paste(
      x$Equation,
      format_number(x$Skewness),
      format_number(x$ExcessKurtosis),
      formatC(x$ShapiroP, format = "e", digits = 2),
      formatC(x$BPP, format = "e", digits = 2),
      format_number(x$FullVarianceLM, 2),
      x$FullVarianceDF,
      formatC(x$FullVarianceP, format = "e", digits = 2),
      sep = " & "
    )
  }, character(1)), " \\\\")
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Baseline-adjusted OLS residual diagnostics. Normality $p$ is from Shapiro--Wilk. Fitted $p$ uses squared residuals regressed on fitted values; full LM is $nR^2$ using all original non-intercept design columns, with the displayed degrees of freedom (df). Tests are diagnostic checks, not tests of causal identification}",
    "\\label{tab:application_diagnostics}",
    "\\small",
    "\\begin{tabular}{lrrrrrrr}",
    "\\toprule",
    "Equation & Skew. & Ex. kurt. & Normality $p$ & Fitted $p$ & Full LM & df & Full $p$ \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/application_diagnostics_table.tex")
}

make_covariate_sensitivity_table <- function() {
  data <- application_covariates[application_covariates$Effect %in% c("PNIE", "TE"), ]
  rows <- vapply(seq_len(nrow(data)), function(i) {
    x <- data[i, ]
    paste0(paste(
      if (x$Specification == "With baseline depression") "Included" else "Excluded",
      method_latex(x$Method), x$Effect, format_number(x$Estimate, 4),
      format_number(x$Lower, 4), format_number(x$Upper, 4), sep = " & "
    ), " \\\\")
  }, character(1))
  write_lines(c(
    "\\begin{table}[htbp]", "\\centering", "\\small",
    "\\caption{JOBS II sensitivity to adjustment for pretreatment depressive symptoms. The sample, remaining covariates and treatment--mediator interaction are unchanged. The specification including baseline depression is primary; all five effects are retained in the accompanying CSV}",
    "\\label{tab:application_covariates}",
    "\\begin{tabular}{lllrrr}", "\\toprule",
    "Baseline depression & Method & Effect & Estimate & Lower & Upper \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"
  ), "results/application_covariate_sensitivity_table.tex")
}

make_application_sensitivity_table <- function() {
  rows <- paste0(vapply(seq_len(nrow(application_sensitivity)), function(i) {
    x <- application_sensitivity[i, ]
    paste(
      format_number(x$StepScale, 2),
      format_number(x$PNIE),
      paste0("[", format_number(x$Lower), ", ", format_number(x$Upper), "]"),
      x$MediatorSelectedStart,
      x$OutcomeSelectedStart,
      x$MediatorCandidateRoots,
      x$OutcomeCandidateRoots,
      formatC(x$MediatorScoreNorm, format = "e", digits = 2),
      formatC(x$OutcomeScoreNorm, format = "e", digits = 2),
      sep = " & "
    )
  }, character(1)), " \\\\")
  write_lines(c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{JOBS II algorithm sensitivity across perturbation scales}",
    "\\label{tab:application_sensitivity}",
    "\\scriptsize",
    "\\begin{tabular}{rrrccccrr}",
    "\\toprule",
    "Scale & PNIE & 95\\% CI & Start M & Start Y & Roots M & Roots Y & Score M & Score Y \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  ), "results/application_sensitivity_table.tex")
}

make_full_table <- function(data, path, caption, label) {
  scenario_order <- c(
    "Gaussian", "Skew-normal", "Asymmetric mixture", "Symmetric bimodal"
  )
  method_order <- c(
    "OLS", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM", "Semiparametric"
  )
  ordered <- data[order(
    data$SampleSize,
    match(data$Scenario, scenario_order),
    match(data$Method, method_order),
    match(data$Effect, c("PNIE", "TNIE", "PNDE", "TNDE", "TE"))
  ), ]
  rows <- paste0(vapply(seq_len(nrow(ordered)), function(i) {
    x <- ordered[i, ]
    paste(
      x$SampleSize,
      scenario_latex_full(x$Scenario),
      method_latex_full(x$Method),
      x$Effect,
      x$Replications,
      format_number(x$Bias),
      format_number(x$RMSE),
      format_number(x$Coverage95),
      format_number(x$MonteCarloSE),
      format_number(x$AvgLength),
      format_number(x$SuccessRate),
      sep = " & "
    )
  }, character(1)), " \\\\")
  if (length(unique(ordered$SampleSize)) > 1L) {
    break_after <- max(which(ordered$SampleSize == ordered$SampleSize[1]))
  } else {
    scenario_index <- match(ordered$Scenario, scenario_order)
    break_after <- max(which(
      scenario_index <= ceiling(length(unique(scenario_index)) / 2)
    ))
  }
  rows <- append(rows, "\\pagebreak", after = break_after)
  write_lines(c(
    "\\begingroup",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{2.1pt}",
    "\\begin{longtable}{rlllrrrrrrr}",
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule",
    "$n$ & Error & Method & Effect & Reps. & Bias & RMSE & Cov. & MCSE & Length & Success \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "$n$ & Error & Method & Effect & Reps. & Bias & RMSE & Cov. & MCSE & Length & Success \\\\",
    "\\midrule",
    "\\endhead",
    rows,
    "\\bottomrule",
    "\\end{longtable}",
    "\\endgroup"
  ), path)
}

make_performance_figure <- function() {
  selected <- main[main$Effect %in% c("PNIE", "TE"), ]
  ols <- selected[selected$Method == "OLS", ]
  proposed <- selected[selected$Method == "Semiparametric", ]
  key <- function(x) paste(x$SampleSize, x$Scenario, x$Effect, sep = "::")
  proposed <- proposed[match(key(ols), key(proposed)), ]
  ratio_rmse <- proposed$RMSE / ols$RMSE
  ratio_length <- proposed$AvgLength / ols$AvgLength
  scenario_order <- c(
    "Gaussian", "Skew-normal", "Asymmetric mixture", "Symmetric bimodal"
  )
  x <- match(ols$Scenario, scenario_order)
  offset <- ifelse(ols$SampleSize == 200, -0.10, 0.10) +
    ifelse(ols$Effect == "PNIE", -0.025, 0.025)
  symbols <- ifelse(ols$Effect == "PNIE", 16, 1)
  line_type <- ifelse(ols$SampleSize == 200, 2, 1)

  grDevices::pdf("figures/simulation_efficiency_ratios.pdf", width = 7.6, height = 3.7)
  old_par <- graphics::par(mfrow = c(1, 2), mar = c(5.2, 4.1, 1.0, 0.5))
  for (values in list(ratio_rmse, ratio_length)) {
    graphics::plot(
      NA,
      xlim = c(0.5, 4.5),
      ylim = range(c(0.2, 1.15, values), finite = TRUE),
      xaxt = "n",
      xlab = "Error distribution",
      ylab = if (identical(values, ratio_rmse)) "RMSE ratio" else "Interval-length ratio",
      bty = "n"
    )
    graphics::axis(1, at = 1:4, labels = c("G", "SN", "AM", "BM"))
    graphics::abline(h = 1, lty = 3, col = "grey55")
    for (i in seq_along(values)) {
      graphics::points(x[i] + offset[i], values[i], pch = symbols[i])
      graphics::segments(
        x[i] + offset[i] - 0.04,
        values[i],
        x[i] + offset[i] + 0.04,
        values[i],
        lty = line_type[i]
      )
    }
  }
  graphics::legend(
    "topright",
    legend = c("PNIE, n=200", "TE, n=200", "PNIE, n=500", "TE, n=500"),
    pch = c(16, 1, 16, 1),
    lty = c(2, 2, 1, 1),
    bty = "n",
    cex = 0.72
  )
  graphics::par(old_par)
  grDevices::dev.off()
}

make_power_figure <- function() {
  grDevices::pdf("figures/size_power_curve.pdf", width = 5.4, height = 3.8)
  old_par <- graphics::par(mar = c(4.2, 4.2, 0.8, 0.6))
  methods <- c("OLS", "Semiparametric")
  graphics::plot(
    NA,
    xlim = range(power$TruePNIE),
    ylim = c(0, 1),
    xlab = "True PNIE",
    ylab = "Empirical rejection probability",
    bty = "n"
  )
  graphics::abline(h = 0.05, lty = 3, col = "grey55")
  for (j in seq_along(methods)) {
    subset <- power[power$Method == methods[j], ]
    subset <- subset[order(subset$TruePNIE), ]
    graphics::lines(
      subset$TruePNIE,
      subset$RejectionRate,
      type = "b",
      pch = c(1, 16)[j],
      lty = c(2, 1)[j]
    )
  }
  graphics::legend(
    "bottomleft",
    legend = c("OLS delta", "Proposed"),
    pch = c(1, 16),
    lty = c(2, 1),
    bty = "n"
  )
  graphics::par(old_par)
  grDevices::dev.off()
}

make_error_figure <- function() {
  set.seed(20260902)
  scenarios <- c(
    "gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture"
  )
  # Generate all panels before opening a device, so a missing dependency cannot
  # silently replace the figure with a partially completed PDF.
  draws <- lapply(scenarios, function(scenario) generate_standardized_error(5000, scenario))
  grDevices::pdf("figures/error_distributions.pdf", width = 7.6, height = 5.2)
  old_par <- graphics::par(
    mfrow = c(2, 2),
    mar = c(2.7, 3.6, 2.0, 0.5),
    oma = c(2.1, 0, 0, 0),
    mgp = c(1.8, 0.55, 0)
  )
  for (j in seq_along(scenarios)) {
    scenario <- scenarios[j]
    draw <- draws[[j]]
    graphics::hist(
      draw,
      breaks = 35,
      probability = TRUE,
      col = "grey85",
      border = "white",
      main = scenario_label(scenario),
      xlab = ""
    )
    graphics::lines(stats::density(draw), lwd = 1.3)
  }
  graphics::mtext("Standardized error", side = 1, outer = TRUE, line = 0.5)
  graphics::par(old_par)
  grDevices::dev.off()
}

make_main_table()
make_comparator_table()
make_application_table()
make_application_full_table()
make_power_table()
make_sensitivity_table()
make_diagnostics_table()
make_regression_table()
make_application_sensitivity_table()
make_covariate_sensitivity_table()
make_full_table(
  main,
  "results/main_simulation_full_table.tex",
  paste(
    "Complete main Monte Carlo results; Reps. is the number of numerically",
    "successful fits used for performance summaries, whereas Success uses all",
    "attempts; errors: G, Gaussian; SN, skew-normal; AM, asymmetric mixture;",
    "BM, symmetric bimodal"
  ),
  "tab:main_simulation_full"
)
make_full_table(
  comparator,
  "results/comparator_full_table.tex",
  paste(
    "Complete comparator Monte Carlo results based on 1,000 attempts per cell; Reps. is the number of numerically",
    "successful fits used for performance summaries, whereas Success uses all",
    "attempts; QB denotes quasi-Bayes and Boot denotes bootstrap; other error",
    "abbreviations are defined in Table S1"
  ),
  "tab:comparator_full"
)
make_performance_figure()
make_power_figure()
make_error_figure()

if (Sys.getenv("JKSS_REUSE_RESULTS", unset = "0") == "1" ||
    Sys.getenv("JKSS_SKIP_SIMULATIONS", unset = "0") == "1") {
  # Reformatting retained results must not relabel their analysis environment.
  writeLines(c("Asset rebuild from retained results; not a new simulation run.",
               utils::capture.output(sessionInfo())),
             "results/asset_sessionInfo.txt")
} else {
  required_jkss_packages(include_comparators = TRUE)
  version_packages <- c("Matrix", "rootSolve", "sn", "mediation", "medflex")
  package_versions <- data.frame(
    Package = version_packages,
    Version = vapply(
      version_packages,
      function(package) as.character(utils::packageVersion(package)),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(package_versions, "results/package_versions.csv", row.names = FALSE)
  writeLines(utils::capture.output(sessionInfo()), "results/sessionInfo.txt")
}
cat("Manuscript assets created; retained-result provenance preserved in reuse mode.\n")

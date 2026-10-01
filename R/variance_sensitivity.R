source(file.path("R", "semiparametric_mediation.R"))

variance_configurations <- function() {
  rbind(
    data.frame(Mechanism = "Constant", Kappa = 0),
    expand.grid(Mechanism = c("Covariate", "Treatment"),
                Kappa = c(0.15, 0.30, 0.60), stringsAsFactors = FALSE)
  )
}

variance_scale <- function(X, T, mechanism, kappa) {
  if (mechanism == "Constant") return(rep(1, length(X)))
  if (mechanism == "Covariate") return(exp(kappa * X / 2 - kappa^2 / 4))
  if (mechanism == "Treatment") {
    return(exp(kappa * (2 * T - 1) / 2) / sqrt(cosh(kappa)))
  }
  stop("Unknown variance mechanism.")
}

variance_design_ratio <- function(mechanism, kappa) {
  if (mechanism == "Covariate") return(exp(diff(qnorm(c(0.25, 0.75))) * kappa))
  if (mechanism == "Treatment") return(exp(2 * kappa))
  1
}

variance_data <- function(X, T, U_M, U_Y, mechanism, kappa) {
  s <- variance_scale(X, T, mechanism, kappa)
  M <- 0.2 + 0.4 * T + 0.3 * X + s * U_M
  Y <- 0.5 * T - 0.8 * M + T * M + 0.4 * X + s * U_Y
  data.frame(Y = Y, M = M, T = T, X = X)
}

variance_task_grid <- function(reps = 500L, seed = 20260917L) {
  grid <- expand.grid(Scenario = c("gaussian", "asymmetric_mixture"),
                      Replication = seq_len(reps), stringsAsFactors = FALSE)
  set.seed(seed)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
  grid$SampleSize <- 300L
  grid
}

variance_worker <- function(task) {
  set.seed(task$Seed)
  n <- task$SampleSize
  X <- rnorm(n)
  T <- rbinom(n, 1, 0.5)
  U_M <- generate_standardized_error(n, task$Scenario)
  U_Y <- generate_standardized_error(n, task$Scenario)
  configs <- variance_configurations()
  records <- status <- list()
  index <- 0L
  for (j in seq_len(nrow(configs))) {
    config <- configs[j, ]
    data <- variance_data(X, T, U_M, U_Y, config$Mechanism, config$Kappa)
    for (method in c("OLS", "Semiparametric")) {
      index <- index + 1L
      failure <- ""
      started <- proc.time()[["elapsed"]]
      fit <- tryCatch(
        fit_stacked_mediation(data, method, "X"),
        error = function(e) {
          failure <<- conditionMessage(e)
          NULL
        }
      )
      valid <- !is.null(fit) && nrow(fit) == 5L &&
        all(is.finite(as.matrix(fit[c("Estimate", "StdError", "Lower", "Upper")]))) &&
        all(fit$StdError > 0)
      if (!valid && !nzchar(failure)) failure <- "Nonfinite or invalid effect output."
      meta <- cbind(task, config, Method = method, Success = valid,
                    FailureReason = failure, stringsAsFactors = FALSE)
      status[[index]] <- cbind(meta, ElapsedSeconds = proc.time()[["elapsed"]] - started)
      if (!valid) {
        fit <- data.frame(Effect = names(true_mediation_effects()), Estimate = NA_real_,
                          StdError = NA_real_, Lower = NA_real_, Upper = NA_real_)
      }
      records[[index]] <- cbind(meta[rep(1, 5), ],
                                fit[c("Effect", "Estimate", "StdError", "Lower", "Upper")])
    }
  }
  list(records = do.call(rbind, records), status = do.call(rbind, status))
}

summarize_variance_records <- function(records) {
  truth <- true_mediation_effects()
  keys <- c("Scenario", "SampleSize", "Mechanism", "Kappa", "Method", "Effect")
  pieces <- split(records, interaction(records[keys], drop = TRUE))
  result <- lapply(pieces, function(d) {
    good <- d[d$Success, ]
    B <- nrow(d)
    N <- nrow(good)
    error <- good$Estimate - truth[d$Effect[1]]
    covered <- good$Lower <= truth[d$Effect[1]] & good$Upper >= truth[d$Effect[1]]
    coverage <- if (N) mean(covered) else NA_real_
    rmse <- if (N) sqrt(mean(error^2)) else NA_real_
    report_cover <- sum(covered) / B
    cbind(d[1, keys], Attempted = B, Valid = N, TrueValue = unname(truth[d$Effect[1]]),
          SuccessRate = N / B, SuccessMCSE = sqrt((N / B) * (1 - N / B) / B),
          Bias = if (N) mean(error) else NA_real_,
          BiasMCSE = if (N > 1) sd(error) / sqrt(N) else NA_real_,
          RMSE = rmse,
          RMSEMCSE = if (N > 1 && rmse > 0) sd(error^2) / sqrt(N) / (2 * rmse) else NA_real_,
          EmpiricalSD = if (N > 1) sd(good$Estimate) else NA_real_,
          MeanSE = if (N) mean(good$StdError) else NA_real_,
          Coverage95 = coverage,
          CoverageMCSE = if (N) sqrt(coverage * (1 - coverage) / N) else NA_real_,
          AvgLength = if (N) mean(good$Upper - good$Lower) else NA_real_,
          ReportAndCover = report_cover,
          ReportAndCoverMCSE = sqrt(report_cover * (1 - report_cover) / B))
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  result[do.call(order, result[keys]), ]
}

paired_variance_summary <- function(records) {
  keys <- c("Scenario", "SampleSize", "Mechanism", "Kappa", "Replication", "Effect")
  paired <- merge(records[records$Method == "OLS" & records$Success, ],
                  records[records$Method == "Semiparametric" & records$Success, ],
                  by = keys, suffixes = c("OLS", "Semi"))
  groups <- c("Scenario", "SampleSize", "Mechanism", "Kappa", "Effect")
  pieces <- split(paired, interaction(paired[groups], drop = TRUE))
  truth <- true_mediation_effects()
  result <- lapply(pieces, function(d) {
    true <- truth[d$Effect[1]]
    mse_ols <- mean((d$EstimateOLS - true)^2)
    mse_semi <- mean((d$EstimateSemi - true)^2)
    cbind(d[1, groups], PairedValid = nrow(d),
          BiasOLS = mean(d$EstimateOLS - true), BiasSemi = mean(d$EstimateSemi - true),
          RMSEOLS = sqrt(mse_ols), RMSESemi = sqrt(mse_semi),
          RMSERatio = sqrt(mse_semi / mse_ols),
          CoverageOLS = mean(d$LowerOLS <= true & d$UpperOLS >= true),
          CoverageSemi = mean(d$LowerSemi <= true & d$UpperSemi >= true),
          LengthRatio = mean(d$UpperSemi - d$LowerSemi) / mean(d$UpperOLS - d$LowerOLS))
  })
  do.call(rbind, result)
}

variance_assets <- function(summary, output_dir = "results", figure_dir = "figures") {
  selected <- summary[summary$Effect %in% c("PNIE", "TE"), ]
  selected <- selected[order(selected$Scenario, selected$Mechanism, selected$Kappa,
                              match(selected$Effect, c("PNIE", "TE")),
                              match(selected$Method, c("OLS", "Semiparametric"))), ]
  lines <- c("\\begingroup\\small\\setlength{\\tabcolsep}{3pt}",
    "\\begin{longtable}{llrllrrrrrr}",
    "\\caption{Variance-misspecification sensitivity at $n=300$, with 500 attempted replications per configuration. Coverage and precision are conditional on numerical success. Bias and coverage MCSEs are in parentheses; Success uses all attempts.}\\label{tab:variance-sensitivity}\\\\",
    "\\toprule", "Error & Scale & $\\kappa$ & Effect & Method & Bias (MCSE) & RMSE & Coverage (MCSE) & Length & Success & R+C\\\\",
    "\\midrule\\endfirsthead", "\\toprule",
    "Error & Scale & $\\kappa$ & Effect & Method & Bias (MCSE) & RMSE & Coverage (MCSE) & Length & Success & R+C\\\\",
    "\\midrule\\endhead", "\\bottomrule\\endfoot")
  previous <- ""
  for (j in seq_len(nrow(selected))) {
    d <- selected[j, ]
    if (d$Scenario != previous && nzchar(previous)) lines <- c(lines, "\\midrule")
    previous <- d$Scenario
    lines <- c(lines, sprintf("%s & %s & %.2f & %s & %s & %.3f (%.3f) & %.3f & %.3f (%.3f) & %.3f & %.3f & %.3f\\\\",
      ifelse(d$Scenario == "gaussian", "G", "AM"),
      switch(d$Mechanism, Constant = "None", Covariate = "$X$", Treatment = "$T$"),
      d$Kappa, d$Effect, ifelse(d$Method == "OLS", "OLS", "Semi"),
      d$Bias, d$BiasMCSE, d$RMSE, d$Coverage95, d$CoverageMCSE, d$AvgLength,
      d$SuccessRate, d$ReportAndCover))
  }
  lines <- c(lines, "\\end{longtable}",
    "\\noindent G: Gaussian; AM: asymmetric mixture; Semi: the unchanged proposed estimator. R+C is the fraction of all attempted data sets yielding an interval covering the truth, not conventional coverage. Scale $X$ denotes covariate-dependent variance and scale $T$ treatment-dependent variance; None is their common constant-variance reference. All five effects and paired-fit summaries are retained in the project repository.",
    "\\endgroup")
  writeLines(lines, file.path(output_dir, "variance_sensitivity_table.tex"))
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(figure_dir, "variance_sensitivity_coverage.pdf"), width = 8, height = 6)
  old <- par(mfrow = c(2, 2), mar = c(3.3, 3.8, 2.2, 0.8), oma = c(4, 0, 0, 0),
             mgp = c(2.2, 0.6, 0), cex = 0.85)
  for (scenario in c("gaussian", "asymmetric_mixture")) {
    for (effect in c("PNIE", "TE")) {
      d <- summary[summary$Scenario == scenario & summary$Effect == effect, ]
      plot(NA, xlim = c(0, 0.6), ylim = c(0, 1), xlab = expression(kappa),
           ylab = "Nominal 95% interval coverage",
           main = paste(ifelse(scenario == "gaussian", "Gaussian", "Asymmetric mixture"), effect))
      abline(h = 0.95, col = "gray55", lty = 3)
      for (mechanism in c("Covariate", "Treatment")) for (method in c("OLS", "Semiparametric")) {
        z <- d[d$Method == method & d$Mechanism %in% c("Constant", mechanism), ]
        z <- z[order(z$Kappa), ]
        lines(z$Kappa, z$Coverage95,
              col = ifelse(mechanism == "Covariate", "#19677B", "#B44E28"),
              lty = ifelse(method == "OLS", 2, 1), pch = ifelse(method == "OLS", 1, 16),
              type = "b", lwd = 1.4)
      }
    }
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom", legend = c("X scale: OLS", "X scale: Semi", "T scale: OLS", "T scale: Semi"),
         col = c("#19677B", "#19677B", "#B44E28", "#B44E28"),
         lty = c(2, 1, 2, 1), pch = c(1, 16, 1, 16), horiz = TRUE, bty = "n", cex = 0.85)
  par(old)
  dev.off()
}

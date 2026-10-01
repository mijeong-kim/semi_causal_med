source(file.path("R", "results_io.R"))
source(file.path("R", "confounding_curve.R"))
out <- Sys.getenv("JKSS_CURVE_OUT", "results")
basis <- read_jkss_csv(file.path(out, "confounding_curve_basis.csv"))
truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.70, TNDE = 1.10, TE = 0.78)
keys <- c("Scenario", "TrueRho", "Method", "Effect")
summaries <- lapply(split(basis, interaction(basis[keys], drop = TRUE)), function(d) {
  good <- d[d$Success, ]
  if (nrow(good) < 2) stop("Fewer than two valid fits in a curve cell.")
  do.call(rbind, lapply(curve_rho_grid(), function(rho) {
    e <- curve_evaluate(good, rho)
    target <- curve_population(d$TrueRho[1], rho)[d$Effect[1]]
    causal <- truth[d$Effect[1]]
    error <- e$Estimate - target
    covered <- e$Lower <= target & e$Upper >= target
    includes_causal <- e$Lower <= causal & e$Upper >= causal
    coverage <- mean(covered)
    cbind(d[1, keys], AssumedRho = rho, CurveTruth = unname(target), CausalTruth = unname(causal),
          Attempted = nrow(d), Valid = nrow(good), SuccessRate = nrow(good) / nrow(d),
          MeanEstimate = mean(e$Estimate), Q05 = unname(quantile(e$Estimate, 0.05)),
          Q95 = unname(quantile(e$Estimate, 0.95)), Bias = mean(error),
          BiasMCSE = sd(error) / sqrt(nrow(good)), RMSE = sqrt(mean(error^2)),
          EmpiricalSD = sd(e$Estimate), MeanSE = mean(e$StdError),
          Coverage95 = coverage, CoverageMCSE = sqrt(coverage * (1 - coverage) / nrow(good)),
          AvgLength = mean(e$Upper - e$Lower), ReportAndCover = sum(covered) / nrow(d),
          CausalBias = mean(e$Estimate - causal),
          CausalInclusion = mean(includes_causal),
          CausalInclusionMCSE = sqrt(mean(includes_causal) * (1 - mean(includes_causal)) / nrow(good)),
          CausalReportAndInclude = sum(includes_causal) / nrow(d))
  }))
})
summary <- do.call(rbind, summaries)
pair_keys <- c("Scenario", "TrueRho", "Replication", "Effect")
pair <- merge(basis[basis$Success & basis$Method == "OLS", ],
              basis[basis$Success & basis$Method == "Semiparametric", ],
              by = pair_keys, suffixes = c("OLS", "Semi"))
paired <- do.call(rbind, lapply(split(pair, interaction(pair[c("Scenario", "TrueRho", "Effect")], drop = TRUE)), function(d) {
  a <- d[paste0(c("A", "B", "VA", "CAB", "VB"), "OLS")]
  b <- d[paste0(c("A", "B", "VA", "CAB", "VB"), "Semi")]
  names(a) <- names(b) <- c("A", "B", "VA", "CAB", "VB")
  do.call(rbind, lapply(curve_rho_grid(), function(rho) {
    x <- curve_evaluate(a, rho)
    y <- curve_evaluate(b, rho)
    target <- curve_population(d$TrueRho[1], rho)[d$Effect[1]]
    sx <- (x$Estimate - target)^2
    sy <- (y$Estimate - target)^2
    ratio <- sqrt(mean(sy) / mean(sx))
    influence <- 0.5 * ratio * ((sy - mean(sy)) / mean(sy) - (sx - mean(sx)) / mean(sx))
    cbind(d[1, c("Scenario", "TrueRho", "Effect")], AssumedRho = rho,
          PairedValid = nrow(d), RMSERatio = ratio, RatioMCSE = sd(influence) / sqrt(nrow(d)),
          CoverageOLS = mean(x$Lower <= target & x$Upper >= target),
          CoverageSemi = mean(y$Lower <= target & y$Upper >= target),
          LengthRatio = mean(y$Upper - y$Lower) / mean(x$Upper - x$Lower))
  }))
}))
write.csv(summary, file.path(out, "confounding_curve_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(out, "confounding_curve_paired.csv"), row.names = FALSE)
thresholds <- read_jkss_csv(file.path(out, "confounding_curve_threshold.csv"))
threshold_summary <- do.call(rbind, lapply(split(thresholds, interaction(thresholds[keys], drop = TRUE)), function(d) {
  good <- d[d$Success, ]
  slope <- if (d$Effect[1] == "PNIE") -0.8 else 0.2
  cv <- slope + d$TrueRho[1]
  target <- cv / sqrt(cv^2 + 1 - d$TrueRho[1]^2)
  error <- good$RhoZero - target
  covered <- good$Lower <= target & good$Upper >= target
  coverage <- mean(covered)
  cbind(d[1, keys], TrueThreshold = target, Attempted = nrow(d), Valid = nrow(good),
        SuccessRate = nrow(good) / nrow(d), Bias = mean(error),
        BiasMCSE = sd(error) / sqrt(nrow(good)), RMSE = sqrt(mean(error^2)),
        Coverage95 = coverage, CoverageMCSE = sqrt(coverage * (1 - coverage) / nrow(good)),
        AvgLength = mean(good$Upper - good$Lower), ReportAndCover = sum(covered) / nrow(d))
}))
write.csv(threshold_summary, file.path(out, "confounding_curve_threshold_summary.csv"), row.names = FALSE)
cat("Full-curve summaries written:", nrow(summary), "cells and", nrow(paired), "paired cells.\n")

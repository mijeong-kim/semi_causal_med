source(file.path("R", "results_io.R"))
out <- Sys.getenv("JKSS_CONFOUNDING_OUT", "results")
records <- read_jkss_csv(file.path(out, "confounding_records.csv"))
thresholds <- read_jkss_csv(file.path(out, "confounding_thresholds.csv"))
truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.7, TNDE = 1.1, TE = 0.78)
metrics <- function(d, parameter, target) {
  good <- d[d$Success, ]
  target <- if (length(target) == 1) rep(target, nrow(good)) else target[d$Success]
  error <- good[[parameter]] - target
  N <- nrow(good); B <- nrow(d)
  coverage <- mean(good$Lower <= target & good$Upper >= target)
  data.frame(Attempted = B, Valid = N, SuccessRate = N / B,
             Bias = mean(error), BiasMCSE = sd(error) / sqrt(N),
             RMSE = sqrt(mean(error^2)), EmpiricalSD = sd(good[[parameter]]),
             MeanSE = mean(good$StdError), Coverage95 = coverage,
             CoverageMCSE = sqrt(coverage * (1 - coverage) / N),
             AvgLength = mean(good$Upper - good$Lower),
             ReportAndCover = sum(good$Lower <= target & good$Upper >= target) / B)
}
keys <- c("Scenario", "TrueRho", "Evaluation", "AssumedRho", "Method", "Effect")
summary <- do.call(rbind, lapply(split(records, interaction(records[keys], drop = TRUE)), function(d) {
  cbind(d[1, keys], metrics(d, "Estimate", truth[d$Effect[1]]))
}))
threshold_keys <- c("Scenario", "TrueRho", "Method", "Effect")
threshold_summary <- do.call(rbind, lapply(split(thresholds, interaction(thresholds[threshold_keys], drop = TRUE)), function(d) {
  cbind(d[1, threshold_keys], metrics(d, "RhoZero", d$TrueThreshold))
}))
pair_keys <- c("Scenario", "TrueRho", "Replication", "Evaluation", "AssumedRho", "Effect")
pair <- merge(records[records$Method == "OLS" & records$Success, ],
              records[records$Method == "Semiparametric" & records$Success, ],
              by = pair_keys, suffixes = c("OLS", "Semi"))
group_keys <- setdiff(pair_keys, "Replication")
paired <- do.call(rbind, lapply(split(pair, interaction(pair[group_keys], drop = TRUE)), function(d) {
  true <- truth[d$Effect[1]]
  cbind(d[1, group_keys], PairedValid = nrow(d),
        RMSERatio = sqrt(mean((d$EstimateSemi - true)^2) / mean((d$EstimateOLS - true)^2)),
        CoverageOLS = mean(d$LowerOLS <= true & d$UpperOLS >= true),
        CoverageSemi = mean(d$LowerSemi <= true & d$UpperSemi >= true),
        LengthRatio = mean(d$UpperSemi - d$LowerSemi) / mean(d$UpperOLS - d$LowerOLS))
}))
for (name in c("summary", "threshold_summary", "paired")) {
  write.csv(get(name), file.path(out, paste0("confounding_", name, ".csv")), row.names = FALSE)
}
cat("Confounding summaries: all effects, two assumed-rho analyses, thresholds and paired fits.\n")

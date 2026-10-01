source(file.path("R", "results_io.R"))
source(file.path("R", "confounding_curve.R"))
out <- Sys.getenv("JKSS_CURVE_OUT", "results")
reps <- as.integer(Sys.getenv("JKSS_CURVE_REPS", "500"))
read_curve <- function(name) read_jkss_csv(file.path(out, paste0("confounding_curve_", name, ".csv")))
basis <- read_curve("basis")
status <- read_curve("status")
threshold <- read_curve("threshold")
summary <- read_curve("summary")
paired <- read_curve("paired")
tsummary <- read_curve("threshold_summary")
seeds <- read_curve("seeds")
id <- c("Scenario", "TrueRho", "Replication", "Method")
stopifnot(nrow(status) == 12 * reps, nrow(basis) == 60 * reps,
          nrow(threshold) == 24 * reps, nrow(summary) == 1980,
          nrow(paired) == 990, nrow(tsummary) == 24,
          !anyDuplicated(status[id]), !anyDuplicated(basis[c(id, "Effect")]),
          all(nzchar(status$FailureReason[!status$Success])),
          all(is.na(basis$A[!basis$Success])),
          all(is.finite(as.matrix(basis[basis$Success, c("A", "B", "VA", "CAB", "VB")]))),
          all(basis$VA[basis$Success] > 0))

# Compare every refit with the original two-value analysis, not only selected examples.
old_status <- read_jkss_csv("results/confounding_status.csv")
matched <- merge(status, old_status, by = id, suffixes = c("New", "Old"))
stopifnot(nrow(matched) == nrow(status), all(matched$SuccessNew == matched$SuccessOld),
          all(matched$FailureReasonNew == matched$FailureReasonOld),
          all(matched$SeedNew == matched$SeedOld))
old_records <- read_jkss_csv("results/confounding_records.csv")
matched <- merge(basis, old_records, by = c(id, "Effect"), suffixes = c("New", "Old"))
good <- matched[matched$SuccessNew, ]
h <- good$AssumedRho / sqrt(1 - good$AssumedRho^2)
estimate <- good$A + h * good$B
se <- sqrt(good$VA + 2 * h * good$CAB + h^2 * good$VB)
refit_error <- max(abs(estimate - good$Estimate), abs(se - good$StdError))
stopifnot(nrow(matched) == 2 * nrow(basis), refit_error < 1e-8)
old_threshold <- read_jkss_csv("results/confounding_thresholds.csv")
matched <- merge(threshold, old_threshold, by = c(id, "Effect"), suffixes = c("New", "Old"))
for (name in c("RhoZero", "StdError", "Lower", "Upper", "MeanMediatorDifference")) {
  stopifnot(max(abs(matched[[paste0(name, "New")]] - matched[[paste0(name, "Old")]]), na.rm = TRUE) < 1e-8)
}
original_seeds <- read_jkss_csv("results/confounding_seeds.csv")
matched <- merge(seeds, original_seeds, by = c("Scenario", "TrueRho", "Replication"))
stopifnot(nrow(matched) == nrow(seeds), all(matched$Seed.x == matched$Seed.y))

truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.70, TNDE = 1.10, TE = 0.78)
keys <- c("Scenario", "TrueRho", "Method", "Effect")
groups <- split(basis, interaction(basis[keys], drop = TRUE))
max_summary_error <- 0
for (i in seq_len(nrow(summary))) {
  z <- summary[i, ]
  key <- as.character(interaction(z[keys], drop = TRUE))
  d <- groups[[key]]
  good <- d[d$Success, ]
  h <- z$AssumedRho / sqrt(1 - z$AssumedRho^2)
  estimate <- good$A + h * good$B
  se <- sqrt(good$VA + 2 * h * good$CAB + h^2 * good$VB)
  delta <- 0.4 * (c(-0.8, 0.2) + z$TrueRho - h * sqrt(1 - z$TrueRho^2))
  target <- setNames(c(delta, 0.78 - delta[2], 0.78 - delta[1], 0.78), names(truth))[z$Effect]
  covered <- abs(estimate - target) <= qnorm(0.975) * se
  causal <- abs(estimate - truth[z$Effect]) <= qnorm(0.975) * se
  expected <- c(CurveTruth = target, MeanEstimate = mean(estimate),
    Bias = mean(estimate - target), RMSE = sqrt(mean((estimate - target)^2)),
    Coverage95 = mean(covered), CoverageMCSE = sqrt(mean(covered) * (1 - mean(covered)) / nrow(good)),
    AvgLength = 2 * qnorm(0.975) * mean(se), ReportAndCover = sum(covered) / reps,
    CausalBias = mean(estimate - truth[z$Effect]), CausalInclusion = mean(causal),
    Q05 = unname(quantile(estimate, 0.05)), Q95 = unname(quantile(estimate, 0.95)))
  names(expected)[1] <- "CurveTruth"
  stopifnot(nrow(d) == reps, z$Valid == nrow(good), z$Attempted == reps)
  max_summary_error <- max(max_summary_error, abs(unlist(z[names(expected)]) - expected))
}
stopifnot(max_summary_error < 1e-9)
for (i in seq_len(nrow(tsummary))) {
  z <- tsummary[i, ]
  keep <- rep(TRUE, nrow(threshold))
  for (key in keys) keep <- keep & threshold[[key]] == z[[key]]
  d <- threshold[keep, ]; good <- d[d$Success, ]
  slope <- if (z$Effect == "PNIE") -0.8 else 0.2
  target <- (slope + z$TrueRho) / sqrt((slope + z$TrueRho)^2 + 1 - z$TrueRho^2)
  covered <- good$Lower <= target & good$Upper >= target
  expected <- c(Bias = mean(good$RhoZero - target),
                RMSE = sqrt(mean((good$RhoZero - target)^2)), Coverage95 = mean(covered),
                AvgLength = mean(good$Upper - good$Lower), ReportAndCover = sum(covered) / reps)
  stopifnot(nrow(d) == reps, nrow(good) == z$Valid,
            max(abs(unlist(z[names(expected)]) - expected)) < 1e-10)
}
correct <- summary[abs(summary$TrueRho - summary$AssumedRho) < 1e-12, ]
stopifnot(max(abs(correct$CurveTruth - correct$CausalTruth)) < 1e-12,
          max(abs(correct$Coverage95 - correct$CausalInclusion)) < 1e-12)

for (i in seq_len(nrow(paired))) {
  z <- paired[i, ]
  a <- groups[[paste(z$Scenario, z$TrueRho, "OLS", z$Effect, sep = ".")]]
  b <- groups[[paste(z$Scenario, z$TrueRho, "Semiparametric", z$Effect, sep = ".")]]
  common <- intersect(a$Replication[a$Success], b$Replication[b$Success])
  a <- a[match(common, a$Replication), ]; b <- b[match(common, b$Replication), ]
  h <- z$AssumedRho / sqrt(1 - z$AssumedRho^2)
  target <- curve_population(z$TrueRho, z$AssumedRho)[z$Effect]
  expected <- sqrt(mean((b$A + h * b$B - target)^2) / mean((a$A + h * a$B - target)^2))
  stopifnot(length(common) == z$PairedValid, abs(expected - z$RMSERatio) < 1e-9)
}
for (d in split(basis[basis$Success, ], interaction(basis[basis$Success, id], drop = TRUE))) {
  for (name in c("A", "B")) {
    v <- setNames(d[[name]], d$Effect)
    stopifnot(abs(v["TE"] - v["PNIE"] - v["TNDE"]) < 1e-10,
              abs(v["TE"] - v["TNIE"] - v["PNDE"]) < 1e-10)
  }
}
report <- c("Full-curve output validation PASSED.",
            paste("Attempted/valid method fits:", nrow(status), sum(status$Success)),
            paste("Maximum original-refit discrepancy:", format(refit_error)),
            paste("Maximum independent curve-summary discrepancy:", format(max_summary_error)),
            "All original seeds, fit statuses, two-value results and zero boundaries reproduced.",
            "Separate curve/causal targets, paired denominators and exact decompositions checked.")
writeLines(report, file.path(out, "confounding_curve_validation.txt"))
cat(paste(report, collapse = "\n"), "\n")

source(file.path("R", "results_io.R"))
source(file.path("R", "variance_sensitivity.R"))
directory <- Sys.getenv("JKSS_VARIANCE_OUT", "results")
expected_reps <- as.integer(Sys.getenv("JKSS_VARIANCE_REPS", "500"))
read_result <- function(stem) read_jkss_csv(file.path(directory, paste0("variance_sensitivity_", stem, ".csv")))
records <- read_result("records")
status <- read_result("status")
summary <- read_result("summary")
paired <- read_result("paired")
seeds <- read_result("seeds")
truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.7, TNDE = 1.1, TE = 0.78)
stopifnot(max(abs(truth - true_mediation_effects())) < 1e-14,
          nrow(seeds) == 2 * expected_reps,
          isTRUE(all.equal(seeds, variance_task_grid(expected_reps), check.attributes = FALSE)),
          nrow(status) == 28 * expected_reps,
          nrow(records) == 140 * expected_reps, nrow(summary) == 140L,
          all(status$SampleSize == 300), all(status$Success[status$Method == "OLS"]),
          !anyDuplicated(status[c("Scenario", "Mechanism", "Kappa", "Replication", "Method")]),
          all(nzchar(status$FailureReason[!status$Success])),
          all(is.na(records$Estimate[!records$Success])),
          all(is.finite(records$Estimate[records$Success])))

for (k in c(0, 0.15, 0.3, 0.6)) {
  # Independently integrate against the known population law, not sample moments.
  ev <- integrate(function(x) exp(k * x - k^2 / 2) * dnorm(x), -10, 10)$value
  et <- mean(exp(k * (2 * c(0, 1) - 1)) / cosh(k))
  stopifnot(abs(ev - 1) < 1e-8, abs(et - 1) < 1e-14)
  x <- c(-1, 0, 1)
  t <- c(0, 1, 0)
  stopifnot(max(abs(variance_scale(x, t, "Covariate", k)^2 -
                     exp(k * x - k^2 / 2))) < 1e-12,
            max(abs(variance_scale(x, t, "Treatment", k)^2 -
                     exp(k * (2 * t - 1)) / cosh(k))) < 1e-12)
}

keys <- c("Scenario", "SampleSize", "Mechanism", "Kappa", "Method", "Effect")
max_metric_error <- 0
for (i in seq_len(nrow(summary))) {
  z <- summary[i, ]
  keep <- rep(TRUE, nrow(records))
  for (key in keys) keep <- keep & records[[key]] == z[[key]]
  d <- records[keep, ]
  good <- d[d$Success, ]
  error <- good$Estimate - truth[z$Effect]
  coverage <- mean(good$Lower <= truth[z$Effect] & good$Upper >= truth[z$Effect])
  expected <- c(Bias = mean(error), RMSE = sqrt(mean(error^2)), Coverage95 = coverage,
                AvgLength = mean(good$Upper - good$Lower),
                SuccessRate = nrow(good) / expected_reps,
                ReportAndCover = sum(good$Lower <= truth[z$Effect] &
                                     good$Upper >= truth[z$Effect]) / expected_reps,
                CoverageMCSE = sqrt(coverage * (1 - coverage) / nrow(good)))
  stopifnot(nrow(d) == expected_reps, nrow(good) == z$Valid)
  max_metric_error <- max(max_metric_error, abs(unlist(z[names(expected)]) - expected))
}
stopifnot(max_metric_error < 1e-10)

good <- records[records$Success, ]
id <- interaction(good[c("Scenario", "Mechanism", "Kappa", "Replication", "Method")], drop = TRUE)
effects <- split(good, id)
decomposition_error <- max(vapply(effects, function(d) {
  value <- setNames(d$Estimate, d$Effect)
  max(abs(c(value["TE"] - value["PNDE"] - value["TNIE"],
            value["TE"] - value["TNDE"] - value["PNIE"])))
}, numeric(1)))
stopifnot(decomposition_error < 1e-10)

for (i in seq_len(nrow(paired))) {
  z <- paired[i, ]
  d <- records[records$Scenario == z$Scenario & records$Mechanism == z$Mechanism &
                 records$Kappa == z$Kappa & records$Effect == z$Effect & records$Success, ]
  semi <- d[d$Method == "Semiparametric", ]
  ols <- d[d$Method == "OLS", ]
  ols <- ols[match(semi$Replication, ols$Replication), ]
  stopifnot(!anyNA(ols$Replication), nrow(semi) == z$PairedValid,
            abs(sqrt(mean((semi$Estimate - truth[z$Effect])^2) /
                       mean((ols$Estimate - truth[z$Effect])^2)) - z$RMSERatio) < 1e-10)
}

# Reconstruct each density's first base data set and check the mean-based map independently.
for (scenario in unique(seeds$Scenario)) {
  task <- seeds[seeds$Scenario == scenario & seeds$Replication == 1L, ]
  set.seed(task$Seed)
  X <- rnorm(300)
  T <- rbinom(300, 1, 0.5)
  U_M <- generate_standardized_error(300, scenario)
  U_Y <- generate_standardized_error(300, scenario)
  configurations <- variance_configurations()
  for (j in seq_len(nrow(configurations))) {
    config <- configurations[j, ]
    d <- variance_data(X, T, U_M, U_Y, config$Mechanism, config$Kappa)
    a <- coef(lm(M ~ T + X, d))
    b <- coef(lm(Y ~ T * M + X, d))
    m0 <- a[1] + a["X"] * mean(d$X)
    independent <- unname(c(a["T"] * b["M"], a["T"] * (b["M"] + b["T:M"]),
      b["T"] + b["T:M"] * m0,
      b["T"] + b["T:M"] * (m0 + a["T"]),
      b["T"] + a["T"] * b["M"] + b["T:M"] * (m0 + a["T"])))
    r <- records[records$Scenario == scenario & records$Replication == 1L &
                   records$Mechanism == config$Mechanism & records$Kappa == config$Kappa &
                   records$Method == "OLS", ]
    r <- r[match(names(truth), r$Effect), ]
    stopifnot(max(abs(independent - r$Estimate)) < 1e-10)
    if (config$Mechanism == "Constant") {
      set.seed(task$Seed)
      original <- generate_mediation_data(300, scenario)
      stopifnot(max(abs(as.matrix(original) - as.matrix(d))) < 1e-12)
    }
  }
}

report <- c("Variance-sensitivity validation PASSED.",
            paste("Attempted method fits:", nrow(status)),
            paste("Valid method fits:", sum(status$Success)),
            paste("Maximum independent summary discrepancy:", format(max_metric_error)),
            paste("Maximum natural-effect decomposition discrepancy:", format(decomposition_error)),
            "Population scale normalization, seed grid, paired denominators and 14 reconstructed OLS maps checked.")
writeLines(report, file.path(directory, "variance_sensitivity_validation.txt"))
cat(paste(report, collapse = "\n"), "\n")

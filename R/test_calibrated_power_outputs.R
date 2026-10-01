source(file.path("R", "results_io.R"))
source(file.path("R", "calibrated_power.R"))
records <- read_jkss_csv("results/calibrated_power_records.csv", stringsAsFactors = FALSE)
seeds <- read_jkss_csv("results/calibrated_power_seeds.csv")
summary <- read_jkss_csv("results/calibrated_power_summary.csv")
critical <- read_jkss_csv("results/calibrated_power_critical.csv")
paired <- read_jkss_csv("results/calibrated_power_paired.csv")
bootstrap <- readRDS("results/calibrated_power_bootstrap.rds")
expected <- power_seed_grid()
close <- function(a, b) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)))
stopifnot(identical(seeds, expected), nrow(records) == 40000L,
          !anyDuplicated(records[c("Task", "Method")]),
          all(table(records$Method) == 10000L),
          all(table(records$Method[records$Stage == "calibration"]) == 5000L),
          all(table(records$Method[records$Stage == "evaluation"],
                    records$Beta2[records$Stage == "evaluation"]) == 1000L),
          nrow(summary) == 20L, nrow(paired) == 15L,
          all(is.na(records$AbsZ[!records$Success])),
          all(is.finite(records$AbsZ[records$Success])),
          all(records$StdError[records$Success] > 0),
          all(is.finite(summary$TwoStageMCSE)),
          all(summary$TwoStageLower <= summary$TwoStageUpper))
close(records$AbsZ[records$Success],
      abs(records$Estimate[records$Success] / records$StdError[records$Success]))
close(records$TruePNIE, -0.26 * records$Beta2)
stopifnot(all(records$ThresholdM[records$Success & startsWith(records$Method, "Huber")] > 0),
          all(records$ThresholdY[records$Success & startsWith(records$Method, "Huber")] > 0))
stopifnot(identical(records$Seed, seeds$Seed[match(records$Task, seeds$Task)]))
for (i in seq_len(nrow(critical))) {
  sample <- records[records$Stage == "calibration" & records$Method == critical$Method[i] & records$Success, ]
  cutoff <- sort(sample$AbsZ)[ceiling(0.95 * (nrow(sample) + 1))]
  close(cutoff, critical$Critical[i])
  close(nrow(sample), critical$Valid[i])
  stopifnot(mean(sample$AbsZ > cutoff) <= 0.05)
}
for (i in seq_len(nrow(summary))) {
  row <- summary[i, ]
  sample <- records[records$Stage == "evaluation" & records$Beta2 == row$Beta2 & records$Method == row$Method, ]
  valid <- sample[sample$Success, ]
  cutoff <- critical$Critical[match(row$Method, critical$Method)]
  close(row$Valid, nrow(valid))
  close(row$CalibratedRate, mean(valid$AbsZ > cutoff))
  close(row$NominalRate, mean(valid$AbsZ > qnorm(0.975)))
  close(row$CalibratedReportReject, sum(valid$AbsZ > cutoff) / nrow(sample))
  close(row$TwoStageMCSE, sd(bootstrap$rate[, i]))
  close(row$TwoStageLower, quantile(bootstrap$rate[, i], 0.025))
  close(row$TwoStageUpper, quantile(bootstrap$rate[, i], 0.975))
  close(row$RMSE, sqrt(mean((valid$Estimate - valid$TruePNIE)^2)))
}
for (i in seq_len(nrow(paired))) {
  row <- paired[i, ]
  sample <- records[records$Stage == "evaluation" & records$Beta2 == row$Beta2, ]
  valid_tasks <- as.integer(names(which(tapply(sample$Success, sample$Task, all))))
  z <- vapply(power_methods, function(m) {
    tmp <- sample[sample$Method == m, ]
    tmp$AbsZ[match(valid_tasks, tmp$Task)]
  }, numeric(length(valid_tasks)))
  cut <- critical$Critical[match(power_methods, critical$Method)]
  rejection <- sweep(z, 2, cut, ">")
  close(row$CommonValid, length(valid_tasks))
  close(row$Difference, mean(rejection[, 4] - rejection[, match(row$Comparator, power_methods)]))
  close(row$TwoStageMCSE, sd(bootstrap$paired[, i]))
}
# Reconstruct data without the study's data-generator wrapper.
audit_tasks <- c(1L, 2500L, 5000L, 5001L, 6000L, 6001L, 7001L, 8001L, 9001L, 10000L)
for (task_id in audit_tasks) {
  task <- seeds[seeds$Task == task_id, ]
  set.seed(task$Seed)
  n <- 300L
  X <- rnorm(n)
  T <- rbinom(n, 1L, 0.5)
  error <- function() {
    p <- c(0.82, 0.18); mu <- c(-0.7, 3.2); sigma <- c(0.35, 0.75)
    component <- sample.int(2L, n, replace = TRUE, prob = p)
    raw <- rnorm(n, mu[component], sigma[component])
    (raw - sum(p * mu)) / sqrt(sum(p * (sigma^2 + (mu - sum(p * mu))^2)))
  }
  em <- error(); ey <- error()
  M <- 0.2 + task$Beta2 * T + 0.3 * X + em
  Y <- 0.5 * T - 0.26 * M + 0.8 * T * M + 0.4 * X + ey
  data <- data.frame(X, T, M, Y)
  a <- lm(M ~ T + X, data); b <- lm(Y ~ T * M + X, data)
  saved <- records[records$Task == task_id & records$Method == "OLS", ]
  close(saved$Estimate, coef(a)["T"] * coef(b)["M"])
  # Also refit all methods to detect method labels, validity and checkpoint drift.
  fresh <- run_power_task(task)
  old <- records[records$Task == task_id, ]
  stopifnot(identical(fresh$Method, old$Method), identical(fresh$Success, old$Success))
  close(fresh$Estimate, old$Estimate)
  close(fresh$StdError, old$StdError)
}
writeLines(c("PASS: 10,000 unique, disjoint-stage task seeds; 40,000 recorded attempts.",
  "PASS: calibration order statistics, held-out summaries and failure denominators.",
  "PASS: paired common-success comparisons and two-stage bootstrap summaries.",
  "PASS: independent data reconstruction and four-method refits on 10 audit tasks.",
  paste("Validation completed:", Sys.time()), capture.output(sessionInfo())),
  "results/calibrated_power_validation.txt")
message("Independent calibrated-power output validation passed.")

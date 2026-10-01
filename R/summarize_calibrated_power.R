source(file.path("R", "results_io.R"))
source(file.path("R", "calibrated_power.R"))
destination <- Sys.getenv("JKSS_POWER_OUTPUT", "results")
records <- read_jkss_csv(file.path(destination, "calibrated_power_records.csv"),
                    stringsAsFactors = FALSE)
grid <- read_jkss_csv(file.path(destination, "calibrated_power_seeds.csv"))
stopifnot(!anyDuplicated(records[c("Task", "Method")]),
          nrow(records) == nrow(grid) * length(power_methods))
statistic_matrix <- function(tasks) {
  result <- matrix(vapply(power_methods, function(method) {
    rows <- records[records$Method == method, ]
    rows$AbsZ[match(tasks, rows$Task)]
  }, numeric(length(tasks))), nrow = length(tasks), ncol = length(power_methods))
  colnames(result) <- power_methods
  result
}
cal <- statistic_matrix(grid$Task[grid$Stage == "calibration"])
beta <- sort(unique(grid$Beta2[grid$Stage == "evaluation"]))
eval <- lapply(beta, function(b) statistic_matrix(grid$Task[grid$Stage == "evaluation" & grid$Beta2 == b]))
critical <- apply(cal, 2L, calibrated_critical)
nominal <- stats::qnorm(0.975)
rate <- function(x, cutoff) mean(x > cutoff, na.rm = TRUE)
report_rate <- function(x, cutoff) mean(is.finite(x) & x > cutoff)
summary <- do.call(rbind, lapply(seq_along(beta), function(i) {
  x <- eval[[i]]
  do.call(rbind, lapply(seq_along(power_methods), function(j) {
    valid <- sum(is.finite(x[, j]))
    pn <- rate(x[, j], nominal)
    pc <- rate(x[, j], critical[j])
    selected <- records$Stage == "evaluation" & records$Beta2 == beta[i] &
      records$Method == power_methods[j] & records$Success
    error <- records$Estimate[selected] - records$TruePNIE[selected]
    data.frame(Beta2 = beta[i], TruePNIE = -0.26 * beta[i], Method = power_methods[j],
      Attempts = nrow(x), Valid = valid, SuccessRate = valid / nrow(x),
      Critical = critical[j], NominalRate = pn, CalibratedRate = pc,
      NominalMCSE = sqrt(pn * (1 - pn) / valid),
      CalibratedMCSE = sqrt(pc * (1 - pc) / valid),
      NominalReportReject = report_rate(x[, j], nominal),
      CalibratedReportReject = report_rate(x[, j], critical[j]),
      Bias = mean(error), RMSE = sqrt(mean(error^2)), row.names = NULL)
  }))
}))
paired <- do.call(rbind, lapply(seq_along(beta), function(i) {
  x <- eval[[i]]
  x <- x[stats::complete.cases(x), , drop = FALSE]
  reject <- sweep(x, 2L, critical, FUN = ">")
  do.call(rbind, lapply(seq_len(3L), function(j) {
    difference <- as.numeric(reject[, 4L]) - as.numeric(reject[, j])
    data.frame(Beta2 = beta[i], Comparator = power_methods[j],
      CommonValid = nrow(x), ProposedRate = mean(reject[, 4L]),
      ComparatorRate = mean(reject[, j]), Difference = mean(difference),
      ConditionalMCSE = stats::sd(difference) / sqrt(length(difference)))
  }))
}))

# Both samples are resampled independently; methods remain paired within a task.
set.seed(20260924)
bootstrap_reps <- 1000L
critical_boot <- matrix(NA_real_, bootstrap_reps, length(power_methods))
rate_boot <- matrix(NA_real_, bootstrap_reps, nrow(summary))
report_boot <- matrix(NA_real_, bootstrap_reps, nrow(summary))
paired_boot <- matrix(NA_real_, bootstrap_reps, nrow(paired))
for (b in seq_len(bootstrap_reps)) {
  index <- sample.int(nrow(cal), nrow(cal), replace = TRUE)
  cutoff <- apply(cal[index, , drop = FALSE], 2L, calibrated_critical)
  critical_boot[b, ] <- cutoff
  for (i in seq_along(beta)) {
    x <- eval[[i]]
    x <- x[sample.int(nrow(x), nrow(x), replace = TRUE), , drop = FALSE]
    columns <- (i - 1L) * 4L + seq_len(4L)
    rate_boot[b, columns] <- vapply(seq_len(4L), function(j) rate(x[, j], cutoff[j]), numeric(1))
    report_boot[b, columns] <- vapply(seq_len(4L), function(j) report_rate(x[, j], cutoff[j]), numeric(1))
    common <- x[stats::complete.cases(x), , drop = FALSE]
    reject <- sweep(common, 2L, cutoff, FUN = ">")
    paired_boot[b, (i - 1L) * 3L + seq_len(3L)] <-
      colMeans(as.numeric(reject[, 4L]) - reject[, seq_len(3L), drop = FALSE])
  }
}
interval <- function(x) apply(x, 2L, stats::quantile, probs = c(0.025, 0.975), names = FALSE)
summary$TwoStageMCSE <- apply(rate_boot, 2L, stats::sd)
summary$TwoStageLower <- interval(rate_boot)[1L, ]
summary$TwoStageUpper <- interval(rate_boot)[2L, ]
summary$ReportTwoStageMCSE <- apply(report_boot, 2L, stats::sd)
summary$ReportTwoStageLower <- interval(report_boot)[1L, ]
summary$ReportTwoStageUpper <- interval(report_boot)[2L, ]
paired$TwoStageMCSE <- apply(paired_boot, 2L, stats::sd)
paired$TwoStageLower <- interval(paired_boot)[1L, ]
paired$TwoStageUpper <- interval(paired_boot)[2L, ]
calibration <- data.frame(Method = power_methods, Attempts = nrow(cal),
  Valid = colSums(is.finite(cal)), Critical = unname(critical),
  CriticalMCSE = apply(critical_boot, 2L, stats::sd),
  CriticalLower = interval(critical_boot)[1L, ], CriticalUpper = interval(critical_boot)[2L, ],
  CalibrationTail = vapply(seq_len(4L), function(j) rate(cal[, j], critical[j]), numeric(1)))
write.csv(summary, file.path(destination, "calibrated_power_summary.csv"), row.names = FALSE)
write.csv(calibration, file.path(destination, "calibrated_power_critical.csv"), row.names = FALSE)
write.csv(paired, file.path(destination, "calibrated_power_paired.csv"), row.names = FALSE)
saveRDS(list(seed = 20260924, replicates = bootstrap_reps, critical = critical_boot,
             rate = rate_boot, report = report_boot, paired = paired_boot,
             method_order = power_methods, beta_order = beta),
        file.path(destination, "calibrated_power_bootstrap.rds"))
signature <- tools::md5sum(c("R/calibrated_power.R", "R/semiparametric_mediation.R",
  "R/run_calibrated_power.R", "POWER_CALIBRATION_PLAN_20260923.md",
  "R/summarize_calibrated_power.R", file.path(destination, "calibrated_power_records.csv")))
writeLines(c(paste("Summarized:", Sys.time()), paste(names(signature), signature),
  "Bootstrap seed: 20260924; replications: 1000; stages independently resampled.",
  capture.output(sessionInfo())), file.path(destination, "calibrated_power_summary_sessionInfo.txt"))
message("Independent calibration, held-out rejection and two-stage uncertainty summarized.")

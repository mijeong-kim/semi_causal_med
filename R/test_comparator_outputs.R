source(file.path("R", "results_io.R"))
# Independently reconstruct the full comparator summaries and representative OLS fits.
source(file.path("R", "semiparametric_mediation.R"))
records <- read_jkss_csv("results/comparator_records.csv")
status <- read_jkss_csv("results/comparator_status.csv")
summary <- read_jkss_csv("results/comparator_summary.csv")
seeds <- read_jkss_csv("results/comparator_seeds.csv")
scenarios <- c("gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture")
methods <- c("OLS", "Semiparametric", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM")
expected_grid <- expand.grid(Scenario = scenarios, Replication = seq_len(1000L),
                            stringsAsFactors = FALSE)
set.seed(20260828)
expected_seeds <- sample.int(.Machine$integer.max, 4000L)
stopifnot(nrow(status) == 20000L, nrow(seeds) == 4000L, nrow(summary) == 100L,
          identical(seeds$Seed, expected_seeds),
          identical(seeds$Scenario, expected_grid$Scenario),
          identical(seeds$Replication, expected_grid$Replication),
          all(seeds$SampleSize == 300L), all(status$SampleSize == 300L),
          setequal(status$Method, methods), !anyNA(status$Success),
          !anyDuplicated(status[c("Scenario", "Replication", "Method")]),
          all(table(status$Scenario, status$Method) == 1000L))
status_key <- function(x) paste(x$Scenario, x$Replication, x$Method, sep = "::")
record_counts <- table(status_key(records))
expected_counts <- ifelse(status$Success, 5L, 0L)
actual_counts <- as.integer(record_counts[status_key(status)])
actual_counts[is.na(actual_counts)] <- 0L
stopifnot(identical(actual_counts, expected_counts),
          all(status_key(records) %in% status_key(status)))
truth <- true_mediation_effects()
max_summary_error <- 0
for (i in seq_len(nrow(summary))) {
  cell <- summary[i, ]
  x <- records[records$Scenario == cell$Scenario & records$Method == cell$Method &
                 records$Effect == cell$Effect, ]
  valid <- status[status$Scenario == cell$Scenario & status$Method == cell$Method, ]
  target <- unname(truth[cell$Effect])
  coverage <- sum(x$Lower <= target & x$Upper >= target) / nrow(x)
  expected <- c(Replications = nrow(x), SuccessRate = sum(valid$Success) / 1000L,
                Bias = mean(x$Estimate) - target,
                RMSE = sqrt(sum((x$Estimate - target)^2) / nrow(x)),
                Coverage95 = coverage, MonteCarloSE = sqrt(coverage * (1 - coverage) / nrow(x)),
                AvgLength = mean(x$Upper - x$Lower))
  max_summary_error <- max(max_summary_error, max(abs(unlist(cell[names(expected)]) - expected)))
}
stopifnot(max_summary_error < 1e-9)
audit_grid <- seeds[seeds$Replication %in% c(1L, 500L, 501L, 1000L), ]
max_ols_error <- 0
for (i in seq_len(nrow(audit_grid))) {
  task <- audit_grid[i, ]
  set.seed(task$Seed)
  data <- generate_mediation_data(300L, task$Scenario)
  fit <- fit_stacked_mediation(data, "OLS", "X")
  old <- records[records$Scenario == scenario_label(task$Scenario) &
                   records$Replication == task$Replication & records$Method == "OLS", ]
  old <- old[match(fit$Effect, old$Effect), ]
  columns <- c("Estimate", "StdError", "Lower", "Upper")
  max_ols_error <- max(max_ols_error, abs(as.matrix(fit[columns]) - as.matrix(old[columns])))
}
stopifnot(is.finite(max_ols_error), max_ols_error < 1e-9)
report <- c(
  "Comparator validation PASSED.",
  "Four error distributions; 1,000 replications per distribution; n=300.",
  paste("Attempted/valid method fits:", nrow(status), sum(status$Success)),
  "All 4,000 task seeds, five-effect availability, and 100 summary cells checked.",
  paste("Maximum independently reconstructed summary discrepancy:", format(max_summary_error, scientific = TRUE)),
  paste("Maximum OLS refit discrepancy on 16 data sets:", format(max_ols_error, scientific = TRUE)),
  paste("Validated:", Sys.time())
)
writeLines(report, "results/comparator_validation.txt")
cat(paste(report, collapse = "\n"), "\n")

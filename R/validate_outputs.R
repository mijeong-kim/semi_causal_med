required_files <- c(
  "main_simulation_records.csv",
  "main_simulation_status.csv",
  "main_simulation_summary.csv",
  "comparator_records.csv",
  "comparator_status.csv",
  "comparator_summary.csv",
  "size_power_summary.csv",
  "algorithm_sensitivity_records.csv",
  "algorithm_sensitivity_summary.csv",
  "application_effects.csv",
  "application_diagnostics.csv"
)
paths <- file.path("results", required_files)
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("Missing output files: ", paste(missing, collapse = ", "))

read_output <- function(name) {
  utils::read.csv(file.path("results", name), check.names = FALSE)
}

main_records <- read_output("main_simulation_records.csv")
main_status <- read_output("main_simulation_status.csv")
main_summary <- read_output("main_simulation_summary.csv")
comparator_records <- read_output("comparator_records.csv")
comparator_status <- read_output("comparator_status.csv")
comparator_summary <- read_output("comparator_summary.csv")
power <- read_output("size_power_summary.csv")
sensitivity_records <- read_output("algorithm_sensitivity_records.csv")
sensitivity <- read_output("algorithm_sensitivity_summary.csv")
application <- read_output("application_effects.csv")
diagnostics <- read_output("application_diagnostics.csv")

expected_effects <- c("PNIE", "TNIE", "PNDE", "TNDE", "TE")
expected_methods <- c(
  "OLS", "Semiparametric", "Imai quasi-Bayes", "Imai bootstrap", "medflex NEM"
)
expected_reps <- c(
  main = 1000L,
  comparator = 500L,
  power = 1000L,
  sensitivity = 300L
)

validate_effect_records <- function(records, key_columns, label) {
  interval_columns <- c("Estimate", "Lower", "Upper")
  if (!all(vapply(records[interval_columns], function(x) all(is.finite(x)), logical(1)))) {
    stop(label, " contains a non-finite estimate or interval endpoint.")
  }
  if (any(records$Lower > records$Upper)) {
    stop(label, " contains an interval with reversed endpoints.")
  }
  key <- do.call(paste, c(records[key_columns], sep = "::"))
  counts <- table(key)
  if (any(counts != length(expected_effects))) {
    stop(label, " has an incomplete five-effect record.")
  }
  order_index <- order(key, match(records$Effect, expected_effects))
  ordered <- records[order_index, ]
  ordered_key <- key[order_index]
  expected_order <- rep(expected_effects, length(unique(ordered_key)))
  if (!identical(as.character(ordered$Effect), expected_order)) {
    stop(label, " has unexpected effect labels or duplicates.")
  }
  estimate <- matrix(
    ordered$Estimate,
    ncol = length(expected_effects),
    byrow = TRUE,
    dimnames = list(NULL, expected_effects)
  )
  max(c(
    abs(estimate[, "TE"] - estimate[, "PNDE"] - estimate[, "TNIE"]),
    abs(estimate[, "TE"] - estimate[, "TNDE"] - estimate[, "PNIE"])
  ))
}

main_decomposition_error <- validate_effect_records(
  main_records,
  c("Study", "SampleSize", "Scenario", "Replication", "Method"),
  "Main simulation"
)
comparator_decomposition_error <- validate_effect_records(
  comparator_records,
  c("Study", "SampleSize", "Scenario", "Replication", "Method"),
  "Comparator simulation"
)
application_decomposition_error <- validate_effect_records(
  application,
  c("Method"),
  "JOBS II application"
)

validate_summary <- function(summary, label) {
  if (!all(summary$Effect %in% expected_effects)) {
    stop(label, " has an unexpected effect label.")
  }
  numeric_columns <- c(
    "Replications", "SuccessRate", "Bias", "RMSE", "Coverage95",
    "MonteCarloSE", "AvgLength"
  )
  if (!all(vapply(summary[numeric_columns], function(x) all(is.finite(x)), logical(1)))) {
    stop(label, " contains non-finite summary values.")
  }
  if (any(summary$SuccessRate < 0 | summary$SuccessRate > 1) ||
      any(summary$Coverage95 < 0 | summary$Coverage95 > 1) ||
      any(summary$RMSE < 0) || any(summary$AvgLength < 0)) {
    stop(label, " contains a metric outside its valid range.")
  }
}

validate_summary(main_summary, "Main summary")
validate_summary(comparator_summary, "Comparator summary")

count_groups <- function(data, columns) {
  length(unique(do.call(paste, c(data[columns], sep = "::"))))
}
main_status_groups <- count_groups(
  main_status,
  c("Study", "SampleSize", "Scenario", "Method")
)
comparator_status_groups <- count_groups(
  comparator_status,
  c("Study", "SampleSize", "Scenario", "Method")
)
if (main_status_groups != 16L) {
  stop("The main status output should contain 16 design-by-method cells.")
}
if (comparator_status_groups != 20L) {
  stop("The comparator status output should contain 20 design-by-method cells.")
}
if (nrow(main_summary) != 5L * main_status_groups) {
  stop("The main summary does not contain five effects for every design cell.")
}
if (nrow(comparator_summary) != 5L * comparator_status_groups) {
  stop("The comparator summary does not contain five effects for every design cell.")
}

main_cell_counts <- table(
  interaction(
    main_status$SampleSize,
    main_status$Scenario,
    main_status$Method,
    drop = TRUE
  )
)
comparator_cell_counts <- table(
  interaction(
    comparator_status$Scenario,
    comparator_status$Method,
    drop = TRUE
  )
)
if (any(main_cell_counts != expected_reps["main"])) {
  stop("At least one main-simulation cell does not contain 1,000 attempts.")
}
if (any(comparator_cell_counts != expected_reps["comparator"])) {
  stop("At least one comparator cell does not contain 500 attempts.")
}
if (!setequal(unique(comparator_status$Method), expected_methods)) {
  stop("The comparator output does not contain all five prespecified methods.")
}

if (!identical(sort(unique(power$Beta2)), c(0, 0.1, 0.2, 0.3, 0.4))) {
  stop("The size-power grid is incomplete.")
}
if (nrow(power) != 10L) stop("The size-power output should contain 10 rows.")
if (!setequal(unique(power$Method), c("OLS", "Semiparametric"))) {
  stop("The size-power output does not contain both prespecified methods.")
}
if (any(power$Replications != expected_reps["power"])) {
  stop("At least one size-power cell does not contain 1,000 attempts.")
}
power_numeric <- c("SuccessRate", "RejectionRate", "MonteCarloSE")
if (!all(vapply(power[power_numeric], function(x) all(is.finite(x)), logical(1)))) {
  stop("The size-power output contains a non-finite summary value.")
}
if (!all(power$SuccessRate >= 0 & power$SuccessRate <= 1) ||
    !all(power$RejectionRate >= 0 & power$RejectionRate <= 1)) {
  stop("The size-power output contains an invalid probability.")
}
if (!identical(sort(unique(sensitivity$StepScale)), c(0.05, 0.1, 0.2))) {
  stop("The algorithm-sensitivity scale grid is incomplete.")
}
if (nrow(sensitivity) != 6L) {
  stop("The algorithm-sensitivity output should contain six rows.")
}
if (length(unique(sensitivity$Scenario)) != 2L) {
  stop("The algorithm-sensitivity output should contain both error scenarios.")
}
if (any(sensitivity$Replications != expected_reps["sensitivity"])) {
  stop("At least one algorithm-sensitivity cell does not contain 300 attempts.")
}
if (nrow(sensitivity_records) != 2L * 3L * expected_reps["sensitivity"]) {
  stop("The algorithm-sensitivity record file should contain 1,800 rows.")
}
sensitivity_cell_counts <- table(
  interaction(
    sensitivity_records$Scenario,
    sensitivity_records$StepScale,
    drop = TRUE
  )
)
if (length(sensitivity_cell_counts) != 6L ||
    any(sensitivity_cell_counts != expected_reps["sensitivity"])) {
  stop("At least one algorithm-sensitivity record cell is incomplete.")
}
sensitivity_numeric <- c(
  "SuccessRate", "MeanPNIE", "SDEstimate",
  "MeanCandidateRootsMediator", "MeanCandidateRootsOutcome"
)
if (!all(vapply(sensitivity[sensitivity_numeric], function(x) all(is.finite(x)), logical(1)))) {
  stop("The algorithm-sensitivity summary contains a non-finite value.")
}
if (any(sensitivity$SuccessRate < 0 | sensitivity$SuccessRate > 1) ||
    any(sensitivity$MeanCandidateRootsMediator < 1) ||
    any(sensitivity$MeanCandidateRootsOutcome < 1)) {
  stop("The algorithm-sensitivity summary contains an invalid metric.")
}
if (!all(c("BPStatistic", "BPP") %in% names(diagnostics))) {
  stop("The application diagnostics do not contain heteroscedasticity checks.")
}
if (!setequal(unique(application$Method), expected_methods)) {
  stop("The application output does not contain all five prespecified methods.")
}

successful_main <- sum(main_status$Success)
successful_comparator <- sum(comparator_status$Success)
if (nrow(main_records) != 5L * successful_main) {
  stop("Main status and record counts disagree.")
}
if (nrow(comparator_records) != 5L * successful_comparator) {
  stop("Comparator status and record counts disagree.")
}

max_decomposition_error <- max(
  main_decomposition_error,
  comparator_decomposition_error,
  application_decomposition_error
)
if (max_decomposition_error >= 1e-8) {
  stop("A total-effect decomposition failed: ", max_decomposition_error)
}

report <- c(
  "Computational Statistics output validation passed.",
  paste("Main status rows:", nrow(main_status)),
  paste("Main successful fits:", successful_main),
  paste("Comparator status rows:", nrow(comparator_status)),
  paste("Comparator successful fits:", successful_comparator),
  paste("Application methods:", length(unique(application$Method))),
  paste("Maximum decomposition error:", format(max_decomposition_error, scientific = TRUE)),
  paste("Validated at:", format(Sys.time(), tz = "Asia/Seoul"))
)
writeLines(report, file.path("results", "validation_report.txt"))
cat(paste(report, collapse = "\n"), "\n")

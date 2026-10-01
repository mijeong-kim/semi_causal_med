source(file.path("R", "results_io.R"))
# Run only the comparator experiment, optionally extending validated retained fits.
source(file.path("R", "run_simulations.R"))
required_jkss_packages(include_comparators = TRUE)
if (!requireNamespace("sn", quietly = TRUE)) stop("Install sn before running the comparator.")
reps <- as_integer_env("JKSS_REPS_COMPARATOR", 1000L)
cores <- as_integer_env("JKSS_CORES", 8L)
stopifnot(is.finite(reps), reps > 0L, is.finite(cores), cores > 0L)
checkpoint_dir <- file.path("tmp", paste0("comparator_checkpoints_", reps))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
reuse_prefix <- Sys.getenv("JKSS_COMPARATOR_REUSE_PREFIX", "0") == "1"
prefix <- NULL
first_rep <- 1L
audit <- character()
seed_grid <- comparator_seed_grid(reps)
source_signature <- tools::md5sum(c("R/run_simulations.R", "R/semiparametric_mediation.R"))
cache_file <- file.path(checkpoint_dir, "validated_prefix.rds")

if (reuse_prefix && file.exists(cache_file)) {
  prefix <- readRDS(cache_file)
  stopifnot(identical(prefix$source_signature, source_signature))
  first_rep <- prefix$reps + 1L
  audit <- prefix$audit
} else if (reuse_prefix) {
  records <- read_jkss_csv("results/comparator_records.csv", stringsAsFactors = FALSE)
  status <- read_jkss_csv("results/comparator_status.csv", stringsAsFactors = FALSE)
  old_reps <- max(status$Replication)
  stopifnot(old_reps < reps, identical(sort(unique(status$Replication)), seq_len(old_reps)),
            all(status$SampleSize == 300L), !anyDuplicated(status[c("Scenario", "Replication", "Method")]),
            all(table(status$Scenario, status$Method) == old_reps))
  old_grid <- comparator_seed_grid(old_reps)
  stopifnot(identical(old_grid$Seed, head(seed_grid$Seed, nrow(old_grid))))
  effect_key <- function(x) paste(x$Scenario, x$Replication, x$Method, x$Effect, sep = "::")
  old_ols <- records[records$Method == "OLS", ]
  ols_lookup <- setNames(old_ols$Estimate, effect_key(old_ols))
  message("Checking all ", nrow(old_grid), " retained data sets against reconstructed OLS effects.")
  discrepancies <- parallel_map(split(old_grid, seq_len(nrow(old_grid))), function(task) {
    set.seed(task$Seed)
    data <- generate_mediation_data(300L, task$Scenario)
    fit <- fit_stacked_mediation(data, "OLS", "X")
    expected <- ols_lookup[paste(scenario_label(task$Scenario), task$Replication,
                                "OLS", fit$Effect, sep = "::")]
    difference <- max(abs(fit$Estimate - expected))
    stopifnot(is.finite(difference), difference < 1e-9)
    difference
  }, cores)
  stopifnot(!any(vapply(discrepancies, inherits, logical(1), "try-error")))
  audit_reps <- unique(c(1L, as.integer(ceiling(old_reps / 2)), old_reps))
  audit_grid <- old_grid[old_grid$Replication %in% audit_reps, , drop = FALSE]
  message("Rechecking all five methods on ", nrow(audit_grid), " retained data sets.")
  checks <- parallel_map(split(audit_grid, seq_len(nrow(audit_grid))), function(task) {
    fresh <- run_comparator_task(task)
    saved_status <- status[status$Scenario == scenario_label(task$Scenario) &
                             status$Replication == task$Replication, ]
    saved_status <- saved_status[match(fresh$status$Method, saved_status$Method), ]
    stopifnot(identical(fresh$status$Success, saved_status$Success))
    saved <- records[match(effect_key(fresh$records), effect_key(records)), ]
    columns <- c("Estimate", "StdError", "Lower", "Upper")
    difference <- max(abs(as.matrix(fresh$records[columns]) - as.matrix(saved[columns])), na.rm = TRUE)
    stopifnot(identical(unname(is.na(fresh$records[columns])), unname(is.na(saved[columns]))),
              is.finite(difference), difference < 1e-9)
    difference
  }, cores)
  stopifnot(!any(vapply(checks, inherits, logical(1), "try-error")))
  audit <- c(
    paste("Retained prefix replications per distribution:", old_reps),
    paste("Reconstructed OLS data sets:", nrow(old_grid)),
    paste("Maximum original OLS discrepancy:", format(max(unlist(discrepancies)), scientific = TRUE)),
    paste("Full five-method audit data sets:", nrow(audit_grid)),
    paste("Maximum five-method discrepancy:", format(max(unlist(checks)), scientific = TRUE))
  )
  prefix <- list(records = records, status = status, reps = old_reps,
                 source_signature = source_signature, audit = audit)
  saveRDS(prefix, cache_file)
  first_rep <- old_reps + 1L
}

message("Comparator: n=300, ", reps, " replications per distribution; quasi-Bayes=1000, bootstrap=499.")
result <- run_comparator_simulation(reps, cores, first_rep = first_rep,
                                   checkpoint_dir = checkpoint_dir)
if (!is.null(prefix)) {
  result$records <- rbind(prefix$records, result$records)
  result$status <- rbind(prefix$status, result$status)
}
stopifnot(nrow(result$status) == 4L * 5L * reps,
          all(table(result$status$Scenario, result$status$Method) == reps),
          !anyDuplicated(result$status[c("Scenario", "Replication", "Method")]))
summary <- summarize_records(result$records, result$status, true_mediation_effects())
dir.create("results", recursive = TRUE, showWarnings = FALSE)
write.csv(result$records, "results/comparator_records.csv", row.names = FALSE)
write.csv(result$status, "results/comparator_status.csv", row.names = FALSE)
write.csv(summary, "results/comparator_summary.csv", row.names = FALSE)
write.csv(seed_grid, "results/comparator_seeds.csv", row.names = FALSE)
writeLines(c(
  paste("Completed:", Sys.time()), "n=300; four error distributions; five methods.",
  paste("Monte Carlo attempts per distribution:", reps),
  "Quasi-Bayesian draws per data set: 1000; percentile-bootstrap resamples: 499.",
  "Seed: 20260828; one stored seed per data set, shared by all five methods.",
  audit,
  paste("New data sets:", 4L * (reps - first_rep + 1L)),
  paste(names(source_signature), source_signature),
  capture.output(sessionInfo())
), "results/comparator_sessionInfo.txt")
message("Comparator results written: ", nrow(result$status), " attempted method fits.")

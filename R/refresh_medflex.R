source(file.path("R", "results_io.R"))
# Repair the natural-effect working model on the identical retained data sets.
# A full future run uses the corrected wrapper directly; this script is a version audit.
source(file.path("R", "run_simulations.R"))
required_jkss_packages(include_comparators = TRUE)
records <- read_jkss_csv("results/comparator_records.csv", stringsAsFactors = FALSE)
status <- read_jkss_csv("results/comparator_status.csv", stringsAsFactors = FALSE)
reps <- max(status$Replication)
grid <- expand.grid(
  Scenario = c("gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture"),
  Replication = seq_len(reps), stringsAsFactors = FALSE
)
set.seed(20260828)
grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
key <- function(x) paste(x$Scenario, x$Replication, x$Effect, sep = "::")
ols_saved <- records[records$Method == "OLS", ]
ols_lookup <- setNames(ols_saved$Estimate, key(ols_saved))
worker <- function(i) {
  task <- grid[i, ]
  set.seed(task$Seed)
  data <- generate_mediation_data(300L, task$Scenario)
  label <- scenario_label(task$Scenario)
  ols <- fit_stacked_mediation(data, "OLS", "X")
  old <- ols_lookup[paste(label, task$Replication, ols$Effect, sep = "::")]
  error <- max(abs(ols$Estimate - old))
  if (!is.finite(error) || error > 1e-9) {
    stop("Stored-seed reconstruction does not reproduce OLS: ", label,
         ", replication ", task$Replication)
  }
  fit <- suppressMessages(fit_medflex(data))
  metadata <- data.frame(Study = "Comparator", SampleSize = 300L,
                         Scenario = label, Replication = task$Replication)
  list(records = effect_records(fit, metadata), ols_error = error)
}
message("Refitting medflex on ", nrow(grid), " seed-reconstructed data sets; checking every OLS match.")
output <- parallel_map(as.list(seq_len(nrow(grid))), worker, cores = 4L)
if (any(vapply(output, inherits, logical(1), "try-error"))) {
  stop("A comparator refresh failed; no retained result files were replaced.")
}
replacement <- do.call(rbind, lapply(output, `[[`, "records"))
stopifnot(nrow(replacement) == 5L * nrow(grid),
          all(is.finite(as.matrix(replacement[c("Estimate", "Lower", "Upper")]))))
index <- which(records$Method == "medflex NEM")
replacement <- replacement[match(key(records[index, ]), key(replacement)), names(records)]
stopifnot(!anyNA(replacement))
records[index, ] <- replacement
status$Success[status$Method == "medflex NEM"] <- TRUE
summary <- summarize_records(records, status, true_mediation_effects())
utils::write.csv(records, "results/comparator_records.csv", row.names = FALSE)
utils::write.csv(status, "results/comparator_status.csv", row.names = FALSE)
utils::write.csv(summary, "results/comparator_summary.csv", row.names = FALSE)
writeLines(c(
  "Only the medflex NEM component was refitted; no new Monte Carlo design or seed was selected.",
  "Imputation: Y ~ T * M + T * X; corrected natural-effect model: Y ~ T0 * T1 + T0 * X.",
  "Covariate X centered at its empirical mean; package-native robust intervals at X=0.",
  paste(nrow(grid), "data sets: four distributions,", reps, "replications each, n=300, seed=20260828."),
  paste("Maximum OLS reconstruction error:",
        format(max(vapply(output, `[[`, numeric(1), "ols_error")), scientific = TRUE)),
  paste("All", nrow(grid), "corrected medflex fits succeeded."),
  utils::capture.output(sessionInfo())
), "results/medflex_refresh_sessionInfo.txt")
message("Corrected medflex records and summaries written; other method records preserved.")

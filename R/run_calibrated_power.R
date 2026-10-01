source(file.path("R", "calibrated_power.R"))
required_jkss_packages()
stopifnot(requireNamespace("quantreg", quietly = TRUE))
integer_env <- function(name, default) {
  value <- as.integer(Sys.getenv(name, as.character(default)))
  stopifnot(length(value) == 1L, is.finite(value), value > 0)
  value
}
calibration <- integer_env("JKSS_POWER_CALIBRATION", 5000L)
evaluation <- integer_env("JKSS_POWER_EVALUATION", 1000L)
cores <- integer_env("JKSS_CORES", 8L)
grid <- power_seed_grid(calibration, evaluation)
destination <- Sys.getenv("JKSS_POWER_OUTPUT", "results")
if ((calibration != 5000L || evaluation != 1000L) && destination == "results") {
  stop("Reduced runs require JKSS_POWER_OUTPUT to protect manuscript results.")
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
checkpoint <- file.path("tmp", paste0("power_checkpoints_", calibration, "_", evaluation))
dir.create(checkpoint, recursive = TRUE, showWarnings = FALSE)
signature <- tools::md5sum(c("R/calibrated_power.R", "R/semiparametric_mediation.R",
                             "R/run_calibrated_power.R", "POWER_CALIBRATION_PLAN_20260923.md"))
metadata <- list(grid = grid, signature = signature,
                 versions = c(R = R.version.string,
                   vapply(c("rootSolve", "Matrix", "quantreg"),
                          function(p) as.character(utils::packageVersion(p)), character(1))))
metadata_file <- file.path(checkpoint, "metadata.rds")
if (file.exists(metadata_file)) stopifnot(identical(readRDS(metadata_file), metadata))
saveRDS(metadata, metadata_file)
write.csv(grid, file.path(destination, "calibrated_power_seeds.csv"), row.names = FALSE)
message("Independent calibration: ", calibration, "; evaluation: ", evaluation,
        " per effect; workers: ", cores)
started <- Sys.time()
files <- file.path(checkpoint, sprintf("task_%05d.rds", grid$Task))
pending <- which(!file.exists(files))
blocks <- split(pending, ceiling(seq_along(pending) / 100L))
for (block in blocks) {
  run_one <- function(i) {
    result <- run_power_task(grid[i, , drop = FALSE])
    temporary <- paste0(files[i], ".part")
    saveRDS(result, temporary)
    if (!file.rename(temporary, files[i])) stop("Could not commit checkpoint.")
    TRUE
  }
  status <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(block, run_one, mc.cores = cores, mc.set.seed = FALSE)
  } else lapply(block, run_one)
  stopifnot(all(vapply(status, identical, logical(1), TRUE)))
  message(sum(file.exists(files)), "/", nrow(grid), " data sets completed; ",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1), " min")
}
records <- do.call(rbind, lapply(files, readRDS))
stopifnot(nrow(records) == length(power_methods) * nrow(grid),
          !anyDuplicated(records[c("Task", "Method")]))
write.csv(records, file.path(destination, "calibrated_power_records.csv"), row.names = FALSE)
writeLines(c(paste("Completed:", Sys.time()), capture.output(str(metadata)),
             capture.output(sessionInfo())),
           file.path(destination, "calibrated_power_sessionInfo.txt"))
message("All independent-calibration experiment records retained.")

source(file.path("R", "results_io.R"))
source(file.path("R", "variance_sensitivity.R"))
reps <- as.integer(Sys.getenv("JKSS_VARIANCE_REPS", "500"))
cores <- as.integer(Sys.getenv("JKSS_VARIANCE_CORES", "4"))
output_dir <- Sys.getenv("JKSS_VARIANCE_OUT", "results")
figure_dir <- if (output_dir == "results") "figures" else file.path(output_dir, "figures")
stopifnot(is.finite(reps), reps > 0, is.finite(cores), cores > 0)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (Sys.getenv("JKSS_VARIANCE_ASSETS_ONLY", "0") == "1") {
  summary <- read_jkss_csv(file.path(output_dir, "variance_sensitivity_summary.csv"))
  variance_assets(summary, output_dir, figure_dir)
  quit(save = "no")
}

tasks <- variance_task_grid(reps)
checkpoint_dir <- file.path("tmp", paste0("variance_checkpoints_", reps))
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
signature <- list(plan = unname(tools::md5sum("VARIANCE_SENSITIVITY_PLAN_20260916.md")),
                  core = unname(tools::md5sum("R/semiparametric_mediation.R")),
                  helpers = unname(tools::md5sum("R/variance_sensitivity.R")), reps = reps)
batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / 24))
all_output <- vector("list", length(batches))
for (b in seq_along(batches)) {
  path <- file.path(checkpoint_dir, sprintf("batch_%03d.rds", b))
  cached <- if (file.exists(path)) readRDS(path) else NULL
  if (!is.null(cached) && identical(cached$signature, signature)) {
    all_output[[b]] <- cached$output
  } else {
    batch <- lapply(batches[[b]], function(i) tasks[i, , drop = FALSE])
    output <- if (.Platform$OS.type == "windows" || cores == 1L) {
      lapply(batch, variance_worker)
    } else {
      parallel::mclapply(batch, variance_worker, mc.cores = cores, mc.set.seed = FALSE)
    }
    if (any(vapply(output, inherits, logical(1), "try-error"))) stop("A worker failed.")
    saveRDS(list(signature = signature, output = output), path)
    all_output[[b]] <- output
  }
  message("Completed base-data tasks ", max(batches[[b]]), "/", nrow(tasks),
          " (seven variance configurations each).")
}
output <- unlist(all_output, recursive = FALSE)
records <- do.call(rbind, lapply(output, `[[`, "records"))
status <- do.call(rbind, lapply(output, `[[`, "status"))
summary <- summarize_variance_records(records)
paired <- paired_variance_summary(records)
for (name in c("records", "status", "summary", "paired")) {
  write.csv(get(name), file.path(output_dir, paste0("variance_sensitivity_", name, ".csv")),
            row.names = FALSE, na = "")
}
write.csv(tasks, file.path(output_dir, "variance_sensitivity_seeds.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "variance_sensitivity_sessionInfo.txt"))
writeLines(capture.output(str(signature)), file.path(output_dir, "variance_sensitivity_source_signature.txt"))
variance_assets(summary, output_dir, figure_dir)
message("Variance sensitivity complete: ", nrow(status), " attempted method fits.")

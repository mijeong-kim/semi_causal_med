source(file.path("R", "confounding_curve.R"))
reps <- as.integer(Sys.getenv("JKSS_CURVE_REPS", "500"))
cores <- as.integer(Sys.getenv("JKSS_CURVE_CORES", "4"))
out <- Sys.getenv("JKSS_CURVE_OUT", "results")
stopifnot(reps > 0, cores > 0)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
grid <- expand.grid(Scenario = c("gaussian", "asymmetric_mixture"), TrueRho = c(-0.3, 0, 0.3),
                    Replication = seq_len(reps), stringsAsFactors = FALSE)
set.seed(20260918)
grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
files <- c("R/confounding_curve.R", "R/run_confounding_curve.R",
           "R/confounding_sensitivity.R", "R/semiparametric_mediation.R",
           "CURVE_SIMULATION_PLAN_20260916.md")
signature <- list(source = tools::md5sum(files), reps = reps, R = R.version.string)
worker <- function(task) {
  set.seed(task$Seed)
  d <- confounding_data(300, task$Scenario, task$TrueRho)
  values <- lapply(c("OLS", "Semiparametric"), function(method) {
    failure <- ""
    result <- tryCatch({
      fit <- fit_reduced_sensitivity(d, method, "X")
      basis <- curve_basis(fit)
      for (rho in curve_rho_grid()) curve_evaluate(basis, rho)
      list(basis = basis, threshold = sensitivity_thresholds(fit))
    }, error = function(e) {
      failure <<- conditionMessage(e)
      NULL
    })
    meta <- cbind(task, Method = method, Success = !is.null(result), FailureReason = failure,
                  stringsAsFactors = FALSE)
    if (is.null(result)) {
      result <- list(basis = data.frame(Effect = names(true_mediation_effects()),
                                        A = NA_real_, B = NA_real_, VA = NA_real_,
                                        CAB = NA_real_, VB = NA_real_),
                     threshold = data.frame(Effect = c("PNIE", "TNIE"), RhoZero = NA_real_,
                       StdError = NA_real_, Lower = NA_real_, Upper = NA_real_,
                       ResidualR2Product = NA_real_, MeanMediatorDifference = NA_real_))
    }
    list(status = meta, basis = cbind(meta[rep(1, 5), ], result$basis),
         threshold = cbind(meta[rep(1, 2), ], result$threshold[setdiff(names(result$threshold), "Method")]))
  })
  lapply(c("status", "basis", "threshold"), function(name) do.call(rbind, lapply(values, `[[`, name)))
}
checkpoint <- file.path("tmp", paste0("curve_checkpoints_", reps))
dir.create(checkpoint, recursive = TRUE, showWarnings = FALSE)
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 24))
all <- vector("list", length(batches))
for (j in seq_along(batches)) {
  path <- file.path(checkpoint, sprintf("batch_%03d.rds", j))
  cache <- if (file.exists(path)) readRDS(path) else NULL
  if (!is.null(cache) && identical(cache$signature, signature)) {
    all[[j]] <- cache$output
  } else {
    tasks <- lapply(batches[[j]], function(i) grid[i, , drop = FALSE])
    output <- if (.Platform$OS.type == "windows" || cores == 1L) lapply(tasks, worker) else
      parallel::mclapply(tasks, worker, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = FALSE)
    if (any(vapply(output, inherits, logical(1), "try-error"))) stop("Curve worker failure.")
    saveRDS(list(signature = signature, output = output), path)
    all[[j]] <- output
  }
  message("Completed curve data sets ", max(batches[[j]]), "/", nrow(grid))
}
all <- unlist(all, recursive = FALSE)
for (k in seq_along(c("status", "basis", "threshold"))) {
  name <- c("status", "basis", "threshold")[k]
  value <- do.call(rbind, lapply(all, `[[`, k))
  write.csv(value, file.path(out, paste0("confounding_curve_", name, ".csv")), row.names = FALSE, na = "")
}
write.csv(grid, file.path(out, "confounding_curve_seeds.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "confounding_curve_sessionInfo.txt"))
writeLines(capture.output(dput(signature)), file.path(out, "confounding_curve_source_signature.txt"))
message("Full-curve refits complete: ", 2 * nrow(grid), " attempted method fits.")

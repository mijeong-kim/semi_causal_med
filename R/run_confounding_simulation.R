source(file.path("R", "confounding_sensitivity.R"))
reps <- as.integer(Sys.getenv("JKSS_CONFOUNDING_REPS", "500"))
cores <- as.integer(Sys.getenv("JKSS_CONFOUNDING_CORES", "4"))
out <- Sys.getenv("JKSS_CONFOUNDING_OUT", "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(reps > 0, cores > 0)
grid <- expand.grid(Scenario = c("gaussian", "asymmetric_mixture"), TrueRho = c(-0.3, 0, 0.3),
                    Replication = seq_len(reps), stringsAsFactors = FALSE)
set.seed(20260918)
grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
signature <- list(source = unname(tools::md5sum(c("R/confounding_sensitivity.R",
                  "R/semiparametric_mediation.R", "CONFOUNDING_SENSITIVITY_PLAN_20260916.md"))),
                  reps = reps, R = R.version.string)
worker <- function(task) {
  set.seed(task$Seed)
  d <- confounding_data(300, task$Scenario, task$TrueRho)
  results <- lapply(c("OLS", "Semiparametric"), function(method) {
    failure <- ""
    fit <- tryCatch(fit_reduced_sensitivity(d, method, "X"), error = function(e) {
      failure <<- conditionMessage(e)
      NULL
    })
    meta <- cbind(task, Method = method, Success = !is.null(fit), FailureReason = failure,
                  stringsAsFactors = FALSE)
    curves <- lapply(c("CorrectRho", "ZeroRho"), function(evaluation) {
      rho <- if (evaluation == "CorrectRho") task$TrueRho else 0
      effects <- if (!is.null(fit)) sensitivity_effects(fit, rho) else
        data.frame(Effect = names(true_mediation_effects()), Estimate = NA_real_,
                   StdError = NA_real_, Lower = NA_real_, Upper = NA_real_, Method = method)
      cbind(meta[rep(1, 5), ], Evaluation = evaluation, AssumedRho = rho,
            effects[setdiff(names(effects), "Method")])
    })
    thresholds <- if (!is.null(fit)) sensitivity_thresholds(fit) else
      data.frame(Effect = c("PNIE", "TNIE"), RhoZero = NA_real_, StdError = NA_real_,
                 Lower = NA_real_, Upper = NA_real_, ResidualR2Product = NA_real_,
                 MeanMediatorDifference = NA_real_, Method = method)
    thresholds$TrueThreshold <- confounding_truth(task$TrueRho)
    list(status = meta, records = do.call(rbind, curves),
         thresholds = cbind(meta[rep(1, 2), ], thresholds[setdiff(names(thresholds), "Method")]))
  })
  lapply(c("status", "records", "thresholds"), function(name) do.call(rbind, lapply(results, `[[`, name)))
}
checkpoints <- file.path("tmp", paste0("confounding_checkpoints_", reps))
dir.create(checkpoints, recursive = TRUE, showWarnings = FALSE)
batches <- split(seq_len(nrow(grid)), ceiling(seq_len(nrow(grid)) / 24))
all <- vector("list", length(batches))
for (b in seq_along(batches)) {
  path <- file.path(checkpoints, sprintf("batch_%03d.rds", b))
  cache <- if (file.exists(path)) readRDS(path) else NULL
  if (!is.null(cache) && identical(cache$signature, signature)) {
    all[[b]] <- cache$output
  } else {
    tasks <- lapply(batches[[b]], function(i) grid[i, , drop = FALSE])
    output <- if (.Platform$OS.type == "windows" || cores == 1) lapply(tasks, worker) else
      parallel::mclapply(tasks, worker, mc.cores = cores, mc.set.seed = FALSE, mc.preschedule = FALSE)
    if (any(vapply(output, inherits, logical(1), "try-error"))) stop("Worker failure.")
    saveRDS(list(signature = signature, output = output), path)
    all[[b]] <- output
  }
  message("Completed confounding data sets ", max(batches[[b]]), "/", nrow(grid))
}
all <- unlist(all, recursive = FALSE)
status <- do.call(rbind, lapply(all, `[[`, 1))
records <- do.call(rbind, lapply(all, `[[`, 2))
thresholds <- do.call(rbind, lapply(all, `[[`, 3))
for (name in c("status", "records", "thresholds")) {
  write.csv(get(name), file.path(out, paste0("confounding_", name, ".csv")), row.names = FALSE, na = "")
}
write.csv(grid, file.path(out, "confounding_seeds.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "confounding_sessionInfo.txt"))
writeLines(capture.output(str(signature)), file.path(out, "confounding_source_signature.txt"))
message("Confounding simulation complete: ", nrow(status), " method attempts.")

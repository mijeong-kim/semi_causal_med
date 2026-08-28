source(file.path("R", "semiparametric_mediation.R"))

as_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else as.integer(value)
}

parallel_map <- function(tasks, fun, cores) {
  if (.Platform$OS.type == "windows" || cores <= 1L) {
    return(lapply(tasks, fun))
  }
  parallel::mclapply(
    tasks,
    fun,
    mc.cores = cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
}

effect_records <- function(table, metadata) {
  if (is.null(table)) {
    return(NULL)
  }
  cbind(metadata, table, stringsAsFactors = FALSE)
}

summarize_records <- function(records, status, truth) {
  status_key <- interaction(
    status$Study,
    status$SampleSize,
    status$Scenario,
    status$Method,
    drop = TRUE
  )
  success_lookup <- tapply(status$Success, status_key, mean)

  record_key <- interaction(
    records$Study,
    records$SampleSize,
    records$Scenario,
    records$Method,
    records$Effect,
    drop = TRUE
  )
  pieces <- split(records, record_key)
  summaries <- lapply(pieces, function(x) {
    true_value <- unname(truth[x$Effect[1]])
    coverage <- mean(x$Lower <= true_value & x$Upper >= true_value)
    status_name <- interaction(
      x$Study[1], x$SampleSize[1], x$Scenario[1], x$Method[1],
      drop = TRUE
    )
    data.frame(
      Study = x$Study[1],
      SampleSize = x$SampleSize[1],
      Scenario = x$Scenario[1],
      Method = x$Method[1],
      Effect = x$Effect[1],
      Replications = nrow(x),
      SuccessRate = unname(success_lookup[as.character(status_name)]),
      Bias = mean(x$Estimate - true_value),
      RMSE = sqrt(mean((x$Estimate - true_value)^2)),
      Coverage95 = coverage,
      MonteCarloSE = sqrt(coverage * (1 - coverage) / nrow(x)),
      AvgLength = mean(x$Upper - x$Lower),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL
  result[order(
    result$Study,
    result$SampleSize,
    result$Scenario,
    result$Method,
    match(result$Effect, names(truth))
  ), ]
}

run_main_simulation <- function(reps, cores, seed = 20260827) {
  sample_sizes <- c(200L, 500L)
  scenarios <- c(
    "gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture"
  )
  grid <- expand.grid(
    SampleSize = sample_sizes,
    Scenario = scenarios,
    Replication = seq_len(reps),
    stringsAsFactors = FALSE
  )
  set.seed(seed)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))

  # Shuffle only the execution order so expensive n = 500 cells are balanced
  # across forked workers; task-specific seeds preserve the simulated data.
  task_order <- sample.int(nrow(grid))
  tasks <- split(grid[task_order, , drop = FALSE], seq_len(nrow(grid)))
  worker <- function(task) {
    set.seed(task$Seed)
    data <- generate_mediation_data(task$SampleSize, task$Scenario)
    method_results <- lapply(c("OLS", "Semiparametric"), function(method) {
      fit <- tryCatch(
        fit_stacked_mediation(data, method, baseline_names = "X"),
        error = function(e) NULL
      )
      metadata <- data.frame(
        Study = "Main",
        SampleSize = task$SampleSize,
        Scenario = scenario_label(task$Scenario),
        Replication = task$Replication,
        stringsAsFactors = FALSE
      )
      list(
        records = effect_records(fit, metadata),
        status = data.frame(
          Study = "Main",
          SampleSize = task$SampleSize,
          Scenario = scenario_label(task$Scenario),
          Replication = task$Replication,
          Method = method,
          Success = !is.null(fit),
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      records = do.call(rbind, lapply(method_results, `[[`, "records")),
      status = do.call(rbind, lapply(method_results, `[[`, "status"))
    )
  }
  output <- parallel_map(tasks, worker, cores)
  list(
    records = do.call(rbind, lapply(output, `[[`, "records")),
    status = do.call(rbind, lapply(output, `[[`, "status"))
  )
}

run_comparator_simulation <- function(reps, cores, seed = 20260828) {
  scenarios <- c(
    "gaussian", "skew_normal", "asymmetric_mixture", "bimodal_mixture"
  )
  grid <- expand.grid(
    Scenario = scenarios,
    Replication = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid$SampleSize <- 300L
  set.seed(seed)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
  tasks <- split(grid, seq_len(nrow(grid)))

  worker <- function(task) {
    set.seed(task$Seed)
    data <- generate_mediation_data(task$SampleSize, task$Scenario)
    fitters <- list(
      "OLS" = function() fit_stacked_mediation(data, "OLS", "X"),
      "Semiparametric" = function() fit_stacked_mediation(data, "Semiparametric", "X"),
      "Imai quasi-Bayes" = function() fit_mediation_package(
        data, boot = FALSE, sims = 1000, seed = task$Seed + 1L
      ),
      "Imai bootstrap" = function() {
        invisible(utils::capture.output(
          value <- suppressMessages(fit_mediation_package(
            data, boot = TRUE, sims = 499, seed = task$Seed + 2L
          ))
        ))
        value
      },
      "medflex NEM" = function() suppressMessages(fit_medflex(data))
    )

    method_results <- lapply(names(fitters), function(method) {
      fit <- tryCatch(
        suppressWarnings(fitters[[method]]()),
        error = function(e) NULL
      )
      metadata <- data.frame(
        Study = "Comparator",
        SampleSize = task$SampleSize,
        Scenario = scenario_label(task$Scenario),
        Replication = task$Replication,
        stringsAsFactors = FALSE
      )
      list(
        records = effect_records(fit, metadata),
        status = data.frame(
          Study = "Comparator",
          SampleSize = task$SampleSize,
          Scenario = scenario_label(task$Scenario),
          Replication = task$Replication,
          Method = method,
          Success = !is.null(fit),
          stringsAsFactors = FALSE
        )
      )
    })
    list(
      records = do.call(rbind, lapply(method_results, `[[`, "records")),
      status = do.call(rbind, lapply(method_results, `[[`, "status"))
    )
  }
  output <- parallel_map(tasks, worker, cores)
  list(
    records = do.call(rbind, lapply(output, `[[`, "records")),
    status = do.call(rbind, lapply(output, `[[`, "status"))
  )
}

run_size_power_simulation <- function(reps, cores, seed = 20260829) {
  beta2_values <- c(0, 0.10, 0.20, 0.30, 0.40)
  gamma <- -0.26
  grid <- expand.grid(
    Beta2 = beta2_values,
    Replication = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid$SampleSize <- 300L
  grid$Scenario <- "asymmetric_mixture"
  set.seed(seed)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
  tasks <- split(grid, seq_len(nrow(grid)))

  worker <- function(task) {
    set.seed(task$Seed)
    data <- generate_mediation_data(
      task$SampleSize,
      task$Scenario,
      beta2 = task$Beta2,
      gamma = gamma,
      eta = 0.8
    )
    method_results <- lapply(c("OLS", "Semiparametric"), function(method) {
      fit <- tryCatch(
        fit_stacked_mediation(data, method, "X"),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        fit <- fit[fit$Effect == "PNIE", ]
      }
      data.frame(
        Beta2 = task$Beta2,
        TruePNIE = task$Beta2 * gamma,
        Replication = task$Replication,
        Method = method,
        Success = !is.null(fit),
        Estimate = if (is.null(fit)) NA_real_ else fit$Estimate,
        Lower = if (is.null(fit)) NA_real_ else fit$Lower,
        Upper = if (is.null(fit)) NA_real_ else fit$Upper,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, method_results)
  }
  records <- do.call(rbind, parallel_map(tasks, worker, cores))
  records$RejectZero <- with(
    records,
    ifelse(Success, Lower > 0 | Upper < 0, NA)
  )
  pieces <- split(records, interaction(records$Beta2, records$Method, drop = TRUE))
  summary <- do.call(rbind, lapply(pieces, function(x) {
    data.frame(
      Beta2 = x$Beta2[1],
      TruePNIE = x$TruePNIE[1],
      Method = x$Method[1],
      Replications = nrow(x),
      SuccessRate = mean(x$Success),
      RejectionRate = mean(x$RejectZero, na.rm = TRUE),
      MonteCarloSE = sqrt(
        mean(x$RejectZero, na.rm = TRUE) *
          (1 - mean(x$RejectZero, na.rm = TRUE)) / sum(x$Success)
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[order(summary$Beta2, summary$Method), ]
}

run_algorithm_sensitivity <- function(reps, cores, seed = 20260830) {
  scales <- c(0.05, 0.10, 0.20)
  grid <- expand.grid(
    Scenario = c("asymmetric_mixture", "bimodal_mixture"),
    Replication = seq_len(reps),
    stringsAsFactors = FALSE
  )
  grid$SampleSize <- 200L
  set.seed(seed)
  grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
  tasks <- split(grid, seq_len(nrow(grid)))

  worker <- function(task) {
    set.seed(task$Seed)
    data <- generate_mediation_data(task$SampleSize, task$Scenario)
    results <- lapply(scales, function(scale) {
      fit <- tryCatch({
        mediator <- fit_semiparametric_lm(M ~ T + X, data, step_scale = scale)
        outcome <- fit_semiparametric_lm(Y ~ T * M + X, data, step_scale = scale)
        stacked <- stack_mediation_fits(mediator, outcome, data["X"])
        effect <- mediation_effects(stacked, "X")
        list(
          effect = effect,
          med_start = mediator$selected_start,
          out_start = outcome$selected_start,
          med_roots = mediator$candidate_count,
          out_roots = outcome$candidate_count
        )
      }, error = function(e) NULL)
      data.frame(
        Scenario = scenario_label(task$Scenario),
        Replication = task$Replication,
        StepScale = scale,
        Success = !is.null(fit),
        PNIE = if (is.null(fit)) NA_real_ else fit$effect$Estimate[fit$effect$Effect == "PNIE"],
        SelectedStartMediator = if (is.null(fit)) NA_integer_ else fit$med_start,
        SelectedStartOutcome = if (is.null(fit)) NA_integer_ else fit$out_start,
        CandidateRootsMediator = if (is.null(fit)) NA_integer_ else fit$med_roots,
        CandidateRootsOutcome = if (is.null(fit)) NA_integer_ else fit$out_roots,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, results)
  }
  records <- do.call(rbind, parallel_map(tasks, worker, cores))
  pieces <- split(records, interaction(records$Scenario, records$StepScale, drop = TRUE))
  summary <- do.call(rbind, lapply(pieces, function(x) {
    data.frame(
      Scenario = x$Scenario[1],
      StepScale = x$StepScale[1],
      Replications = nrow(x),
      SuccessRate = mean(x$Success),
      MeanPNIE = mean(x$PNIE, na.rm = TRUE),
      SDEstimate = stats::sd(x$PNIE, na.rm = TRUE),
      MeanCandidateRootsMediator = mean(x$CandidateRootsMediator, na.rm = TRUE),
      MeanCandidateRootsOutcome = mean(x$CandidateRootsOutcome, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  list(records = records, summary = summary)
}

reps_main <- as_integer_env("CS_REPS_MAIN", 1000L)
reps_comparator <- as_integer_env("CS_REPS_COMPARATOR", 500L)
reps_power <- as_integer_env("CS_REPS_POWER", 1000L)
reps_sensitivity <- as_integer_env("CS_REPS_SENSITIVITY", 300L)
detected_cores <- parallel::detectCores(logical = FALSE)
if (!is.finite(detected_cores)) {
  detected_cores <- parallel::detectCores(logical = TRUE)
}
if (!is.finite(detected_cores)) {
  detected_cores <- 2L
}
default_cores <- max(1L, min(8L, as.integer(detected_cores) - 1L))
cores <- as_integer_env("CS_CORES", default_cores)
if (!is.finite(cores) || cores < 1L) {
  stop("CS_CORES must be a positive integer.")
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
cat("Running main simulation with", reps_main, "replications and", cores, "cores\n")
main <- run_main_simulation(reps_main, cores)
truth <- true_mediation_effects()
main_summary <- summarize_records(main$records, main$status, truth)
utils::write.csv(main$records, "results/main_simulation_records.csv", row.names = FALSE)
utils::write.csv(main$status, "results/main_simulation_status.csv", row.names = FALSE)
utils::write.csv(main_summary, "results/main_simulation_summary.csv", row.names = FALSE)
cat("Main simulation checkpoint written to results/.\n")

cat("Running comparator simulation with", reps_comparator, "replications\n")
comparator <- run_comparator_simulation(reps_comparator, cores)
comparator_summary <- summarize_records(comparator$records, comparator$status, truth)
utils::write.csv(comparator$records, "results/comparator_records.csv", row.names = FALSE)
utils::write.csv(comparator$status, "results/comparator_status.csv", row.names = FALSE)
utils::write.csv(comparator_summary, "results/comparator_summary.csv", row.names = FALSE)
cat("Comparator simulation checkpoint written to results/.\n")

cat("Running size and power simulation with", reps_power, "replications\n")
power <- run_size_power_simulation(reps_power, cores)
utils::write.csv(power, "results/size_power_summary.csv", row.names = FALSE)
cat("Size and power checkpoint written to results/.\n")

cat("Running algorithm sensitivity study with", reps_sensitivity, "replications\n")
sensitivity <- run_algorithm_sensitivity(reps_sensitivity, cores)
utils::write.csv(sensitivity$records, "results/algorithm_sensitivity_records.csv", row.names = FALSE)
utils::write.csv(sensitivity$summary, "results/algorithm_sensitivity_summary.csv", row.names = FALSE)
cat("Algorithm sensitivity checkpoint written to results/.\n")
cat("All simulation outputs written to results/.\n")

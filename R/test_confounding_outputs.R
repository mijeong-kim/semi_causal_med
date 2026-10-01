source(file.path("R", "results_io.R"))
source(file.path("R", "confounding_sensitivity.R"))
directory <- Sys.getenv("JKSS_CONFOUNDING_OUT", "results")
reps <- as.integer(Sys.getenv("JKSS_CONFOUNDING_REPS", "500"))
read_result <- function(stem) read_jkss_csv(file.path(directory, paste0("confounding_", stem, ".csv")))
records <- read_result("records")
status <- read_result("status")
thresholds <- read_result("thresholds")
summary <- read_result("summary")
threshold_summary <- read_result("threshold_summary")
paired <- read_result("paired")
seeds <- read_result("seeds")
truth <- c(PNIE = -0.32, TNIE = 0.08, PNDE = 0.7, TNDE = 1.1, TE = 0.78)
grid <- expand.grid(Scenario = c("gaussian", "asymmetric_mixture"), TrueRho = c(-0.3, 0, 0.3),
                    Replication = seq_len(reps), stringsAsFactors = FALSE)
set.seed(20260918)
grid$Seed <- sample.int(.Machine$integer.max, nrow(grid))
id <- c("Scenario", "TrueRho", "Replication", "Method")
stopifnot(isTRUE(all.equal(grid, seeds, check.attributes = FALSE)),
          nrow(status) == 12 * reps, nrow(records) == 120 * reps,
          nrow(thresholds) == 24 * reps, nrow(summary) == 120,
          nrow(threshold_summary) == 24, nrow(paired) == 60,
          !anyDuplicated(status[id]),
          !anyDuplicated(records[c(id, "Evaluation", "Effect")]),
          !anyDuplicated(thresholds[c(id, "Effect")]),
          all(nzchar(status$FailureReason[!status$Success])),
          all(is.na(records$Estimate[!records$Success])),
          all(is.finite(records$Estimate[records$Success])),
          all(is.finite(records$StdError[records$Success])),
          all(records$StdError[records$Success] > 0))

select_cell <- function(d, z, keys) {
  keep <- rep(TRUE, nrow(d))
  for (key in keys) keep <- keep & d[[key]] == z[[key]]
  d[keep, , drop = FALSE]
}
max_metric_error <- 0
for (kind in c("effects", "thresholds")) {
  tab <- if (kind == "effects") summary else threshold_summary
  raw <- if (kind == "effects") records else thresholds
  keys <- c("Scenario", "TrueRho", "Method", "Effect")
  if (kind == "effects") keys <- c(keys, "Evaluation", "AssumedRho")
  for (i in seq_len(nrow(tab))) {
    z <- tab[i, ]
    d <- select_cell(raw, z, keys)
    good <- d[d$Success, ]
    target <- if (kind == "effects") truth[z$Effect] else good$TrueThreshold
    estimate <- if (kind == "effects") good$Estimate else good$RhoZero
    error <- estimate - target
    covered <- good$Lower <= target & good$Upper >= target
    coverage <- mean(covered)
    expected <- c(Bias = mean(error), BiasMCSE = sd(error) / sqrt(nrow(good)),
                  RMSE = sqrt(mean(error^2)), EmpiricalSD = sd(estimate),
                  MeanSE = mean(good$StdError), Coverage95 = coverage,
                  CoverageMCSE = sqrt(coverage * (1 - coverage) / nrow(good)),
                  AvgLength = mean(good$Upper - good$Lower),
                  SuccessRate = nrow(good) / reps, ReportAndCover = sum(covered) / reps)
    stopifnot(nrow(d) == reps, z$Attempted == reps, z$Valid == nrow(good))
    max_metric_error <- max(max_metric_error, abs(unlist(z[names(expected)]) - expected))
  }
}
stopifnot(max_metric_error < 1e-10)

good <- records[records$Success, ]
groups <- split(good, interaction(good[c(id, "Evaluation")], drop = TRUE))
decomposition_error <- max(vapply(groups, function(d) {
  b <- setNames(d$Estimate, d$Effect)
  max(abs(c(b["TE"] - b["PNDE"] - b["TNIE"], b["TE"] - b["TNDE"] - b["PNIE"])))
}, numeric(1)))
stopifnot(decomposition_error < 1e-10)
zero_true <- records[records$TrueRho == 0, ]
same <- merge(zero_true[zero_true$Evaluation == "CorrectRho", ],
              zero_true[zero_true$Evaluation == "ZeroRho", ], by = c(id, "Effect"))
stopifnot(isTRUE(all.equal(same$Estimate.x, same$Estimate.y)))

for (i in seq_len(nrow(paired))) {
  z <- paired[i, ]
  d <- select_cell(good, z, c("Scenario", "TrueRho", "Evaluation", "AssumedRho", "Effect"))
  a <- d[d$Method == "OLS", ]
  b <- d[d$Method == "Semiparametric", ]
  common <- intersect(a$Replication, b$Replication)
  a <- a[match(common, a$Replication), ]
  b <- b[match(common, b$Replication), ]
  expected <- c(RMSERatio = sqrt(mean((b$Estimate - truth[z$Effect])^2) /
                                  mean((a$Estimate - truth[z$Effect])^2)),
                LengthRatio = mean(b$Upper - b$Lower) / mean(a$Upper - a$Lower),
                CoverageOLS = mean(a$Lower <= truth[z$Effect] & a$Upper >= truth[z$Effect]),
                CoverageSemi = mean(b$Lower <= truth[z$Effect] & b$Upper >= truth[z$Effect]))
  stopifnot(length(common) == z$PairedValid,
            max(abs(unlist(z[names(expected)]) - expected)) < 1e-10)
}

# Reconstruct six data sets, using lm and direct moment algebra rather than the fitted map.
for (i in which(seeds$Replication == 1)) {
  z <- seeds[i, ]
  set.seed(z$Seed)
  d <- confounding_data(300, z$Scenario, z$TrueRho)
  m <- lm(M ~ T + X, d)
  w <- lapply(0:1, function(t) lm(Y ~ X, d[d$T == t, ]))
  vm <- mean(residuals(m)^2)
  vw <- vapply(w, function(f) mean(residuals(f)^2), numeric(1))
  cv <- vapply(0:1, function(t) mean(residuals(m)[d$T == t] * residuals(w[[t + 1]])), numeric(1))
  total <- sum((coef(w[[2]]) - coef(w[[1]])) * c(1, mean(d$X)))
  for (evaluation in c("CorrectRho", "ZeroRho")) {
    rho <- if (evaluation == "CorrectRho") z$TrueRho else 0
    delta <- unname(coef(m)["T"]) * (cv / vm - rho / sqrt(1 - rho^2) * sqrt(vw / vm - cv^2 / vm^2))
    expected <- c(delta, total - delta[2], total - delta[1], total)
    r <- select_cell(records, z, c("Scenario", "TrueRho", "Replication"))
    r <- r[r$Method == "OLS" & r$Evaluation == evaluation, ]
    r <- r[match(names(truth), r$Effect), ]
    stopifnot(max(abs(expected - r$Estimate)) < 1e-10)
  }
  r <- select_cell(thresholds, z, c("Scenario", "TrueRho", "Replication"))
  r <- r[r$Method == "OLS", ]
  r <- r[match(c("PNIE", "TNIE"), r$Effect), ]
  pop_cv <- c(-0.8, 0.2) + z$TrueRho
  stopifnot(max(abs(cv / sqrt(vm * vw) - r$RhoZero)) < 1e-10,
            max(abs(pop_cv / sqrt(pop_cv^2 + 1 - z$TrueRho^2) - r$TrueThreshold)) < 1e-10)
}

if (directory == "results") {
  curves <- read_result("application_curves")
  contour <- read_result("application_contours")
  boundary <- read_result("application_thresholds")
  stopifnot(nrow(curves) == 810, nrow(contour) == 6724, nrow(boundary) == 4,
            max(abs(contour$Rho - contour$Sign * sqrt(contour$R2M * contour$R2Y))) < 1e-12,
            max(abs(boundary$ResidualR2Product - boundary$RhoZero^2)) < 1e-12,
            all(boundary$Lower > -1 & boundary$Upper < 1),
            all(is.finite(curves$StdError)), all(curves$StdError > 0))
  for (d in split(curves, interaction(curves[c("Method", "Rho")], drop = TRUE))) {
    e <- setNames(d$Estimate, d$Effect)
    stopifnot(abs(e["TE"] - e["PNIE"] - e["TNDE"]) < 1e-10,
              abs(e["TE"] - e["TNIE"] - e["PNDE"]) < 1e-10)
  }
  for (method in unique(curves$Method)) {
    for (effect in c("PNIE", "TNIE")) {
      d <- curves[curves$Method == method & curves$Effect == effect, ]
      d <- d[order(d$Rho), ]
      root <- approx(d$Estimate, d$Rho, xout = 0)$y
      reported <- boundary$RhoZero[boundary$Method == method & boundary$Effect == effect]
      stopifnot(abs(root - reported) < 0.001)
    }
    d <- curves[curves$Method == method & curves$Effect == "TE", ]
    stopifnot(diff(range(d$Estimate)) < 1e-12, diff(range(d$StdError)) < 1e-12)
  }
}
report <- c("Confounding-sensitivity output validation PASSED.",
            paste("Attempted/valid method fits:", nrow(status), sum(status$Success)),
            paste("Maximum independent summary discrepancy:", format(max_metric_error)),
            paste("Maximum decomposition discrepancy:", format(decomposition_error)),
            "Seed grid, paired denominators, six reconstructed OLS fits, population zero boundaries checked.",
            if (directory == "results") "Application curves, both contour signs and zero boundaries checked.")
writeLines(report, file.path(directory, "confounding_validation.txt"))
cat(paste(report, collapse = "\n"), "\n")

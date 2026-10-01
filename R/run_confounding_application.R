source(file.path("R", "confounding_sensitivity.R"))
data_environment <- new.env(parent = emptyenv())
load(file.path("data", "jobs.RData"), envir = data_environment)
jobs <- data_environment$jobs
data <- data.frame(
  Y = jobs$depress2,
  M = jobs$job_seek,
  T = jobs$treat,
  depress1 = jobs$depress1,
  econ_hard = jobs$econ_hard,
  sex = jobs$sex,
  age = jobs$age
)
baseline <- c("depress1", "econ_hard", "sex", "age")
stopifnot(nrow(data) == 899L, !anyNA(data))
fits <- lapply(c("OLS", "Semiparametric"), function(method) fit_reduced_sensitivity(data, method, baseline))
names(fits) <- c("OLS", "Semiparametric")
curves <- do.call(rbind, lapply(fits, function(fit) do.call(rbind, lapply(seq(-0.8, 0.8, by = 0.02), function(rho) {
  cbind(Rho = rho, sensitivity_effects(fit, rho))
}))))
thresholds <- do.call(rbind, lapply(fits, sensitivity_thresholds))
grid <- expand.grid(R2M = seq(0, 0.8, length.out = 41), R2Y = seq(0, 0.8, length.out = 41))
contours <- do.call(rbind, lapply(fits, function(fit) do.call(rbind, lapply(c(-1, 1), function(sign) {
  values <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
    rho <- sign * sqrt(grid$R2M[j] * grid$R2Y[j])
    effect <- sensitivity_effects(fit, rho)
    cbind(grid[j, ], Sign = sign, Rho = rho, effect[effect$Effect == "PNIE", ])
  }))
  values
}))))
diagnostics <- do.call(rbind, lapply(fits, function(fit) do.call(rbind, lapply(names(fit$fits), function(name) {
  f <- fit$fits[[name]]
  e <- f$response - drop(f$design %*% f$coefficients)
  auxiliary <- lm(e^2 ~ f$design[, -1, drop = FALSE])
  statistic <- length(e) * summary(auxiliary)$r.squared
  df <- ncol(f$design) - 1
  data.frame(Method = fit$method, Equation = name, N = length(e),
             VarianceLM = statistic, DF = df, ReferenceP = pchisq(statistic, df, lower.tail = FALSE),
             Skewness = mean(e^3) / mean(e^2)^1.5,
             ExcessKurtosis = mean(e^4) / mean(e^2)^2 - 3,
             ScoreNorm = max(abs(colMeans(f$score))))
}))))
for (name in c("curves", "thresholds", "contours", "diagnostics")) {
  write.csv(get(name), file.path("results", paste0("confounding_application_", name, ".csv")), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), "results/confounding_application_sessionInfo.txt")
writeLines(c("JOBS II: same 899 observations and baseline terms as the primary analysis.",
             "Pooled M ~ T + baseline; separate reduced Y ~ baseline in T=0 and T=1.",
             "The reduced-form extension estimates Y models separately by treatment stratum.",
             "A common analyst-specified rho is varied over [-0.8,0.8].",
             "Contour axes are residual R-squared fractions under the shared-confounder calibration.",
             "Both confounding signs and pointwise interval-zero contours are shown.",
             "The reduced-form results complement the primary conditional-outcome analysis.",
             "Working-model variance restrictions and causal qualifications apply."),
           "results/confounding_application_specification.txt")
cat("Confounding application complete: 899 participants, two methods, both confounding signs.\n")

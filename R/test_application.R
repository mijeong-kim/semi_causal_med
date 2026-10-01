source(file.path("R", "results_io.R"))
data_environment <- new.env(parent = emptyenv())
load(file.path("data", "jobs.RData"), envir = data_environment)
jobs <- data_environment$jobs
data <- data.frame(M = jobs$job_seek, Y = jobs$depress2, T = jobs$treat,
                   depress1 = jobs$depress1, econ_hard = jobs$econ_hard,
                   sex = jobs$sex, age = jobs$age)
stopifnot(nrow(data) == 899L, !anyNA(data),
          identical(as.integer(table(data$T)), c(299L, 600L)))
fits <- list(
  Mediator = lm(M ~ T + depress1 + econ_hard + sex + age, data = data),
  Outcome = lm(Y ~ T * M + depress1 + econ_hard + sex + age, data = data)
)
diagnostics <- read_jkss_csv("results/application_diagnostics.csv")
coefficients <- read_jkss_csv("results/application_regression.csv")
for (equation in names(fits)) {
  fit <- fits[[equation]]
  saved <- diagnostics[diagnostics$Equation == equation, ]
  b <- coefficients[coefficients$Equation == equation & coefficients$Method == "OLS", ]
  stopifnot(max(abs(b$Estimate - coef(fit)[b$Term])) < 1e-9)
  x <- model.matrix(fit)
  # Independent QR projection of centered squared residuals checks the LM calculation.
  z <- residuals(fit)^2 - mean(residuals(fit)^2)
  q <- qr.Q(qr(x))
  statistic <- length(z) * sum(crossprod(q, z)^2) / sum(z^2)
  stopifnot(abs(statistic - saved$FullVarianceLM) < 1e-9,
            saved$FullVarianceDF == ncol(x) - 1L,
            abs(saved$FullVarianceP - pchisq(statistic, ncol(x) - 1L,
                                           lower.tail = FALSE)) < 1e-12)
}
support <- read_jkss_csv("results/application_support.csv")
for (i in seq_len(nrow(support))) {
  x <- data[[support$Variable[i]]]
  stopifnot(support$UniqueValues[i] == length(unique(x)),
            abs(support$AtMinimum[i] - mean(x == min(x))) < 1e-12,
            abs(support$AtMaximum[i] - mean(x == max(x))) < 1e-12)
}
effects <- read_jkss_csv("results/application_effects.csv")
sensitivity <- read_jkss_csv("results/application_covariate_sensitivity.csv")
main <- sensitivity[sensitivity$Specification == "With baseline depression", ]
key <- function(x) paste(x$Method, x$Effect, sep = "::")
index <- match(key(main), key(effects))
stopifnot(!anyNA(index), nrow(sensitivity) == 20L,
          max(abs(main$Estimate - effects$Estimate[index])) < 1e-12,
          max(abs(main$Lower - effects$Lower[index])) < 1e-12,
          max(abs(main$Upper - effects$Upper[index])) < 1e-12,
          sum(coefficients$Term == "depress1") == 4L,
          all(effects$SampleSize == nrow(data)))
cat("Application checks passed: JOBS II source, baseline adjustment, diagnostics and sensitivity.\n")

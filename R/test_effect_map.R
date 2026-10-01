source(file.path("R", "semiparametric_mediation.R"))
set.seed(20260915)
n <- 180L
X <- rnorm(n, mean = 0.8)
T <- rbinom(n, 1, 0.5)
M <- 0.2 + 0.4 * T + 0.3 * X + rnorm(n)
Y <- 0.5 * T - 0.8 * M + T * M + 0.4 * X + rnorm(n)
d <- data.frame(X, T, M, Y)
mf <- fit_ols_lm(M ~ T + X, d)
yf <- fit_ols_lm(Y ~ T * M + X, d)
fit <- stack_mediation_fits(mf, yf, d["X"])
out <- mediation_effects(fit, "X")
effect_map <- function(theta) {
  mu0 <- theta["m:(Intercept)"] + theta["m:X"] * theta["xbar:X"]
  a <- theta["m:T"]; b <- theta["y:M"]; interaction <- theta["y:T:M"]
  direct <- theta["y:T"]
  unname(c(a*b, a*(b+interaction), direct+interaction*mu0,
           direct+interaction*(mu0+a), direct+a*b+interaction*(mu0+a)))
}
theta <- fit$parameter
stopifnot(max(abs(out$Estimate - effect_map(theta))) < 1e-12)
G <- sapply(seq_along(theta), function(j) {
  h <- 1e-5 * max(1, abs(theta[j])); step <- theta * 0; step[j] <- h
  (effect_map(theta + step) - effect_map(theta - step)) / (2*h)
})
V <- G %*% fit$covariance %*% t(G)
stopifnot(max(abs(out$StdError^2 - diag(V))) < 1e-9)
constraints <- rbind(c(0, -1, -1, 0, 1), c(-1, 0, 0, -1, 1))
stopifnot(max(abs(constraints %*% G)) < 1e-9,
          max(abs(constraints %*% V)) < 1e-9)
cat("Effect map, baseline-mean propagation and both covariance constraints passed.\n")

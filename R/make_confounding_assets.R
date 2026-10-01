source(file.path("R", "results_io.R"))
summary <- read_jkss_csv("results/confounding_summary.csv")
selected <- summary[summary$Evaluation == "CorrectRho" & summary$Effect %in% c("PNIE", "TE"), ]
selected <- selected[order(selected$Scenario, selected$TrueRho, selected$Effect, selected$Method), ]
lines <- c("\\begin{table}[htbp]\\centering\\small\\setlength{\\tabcolsep}{3pt}",
  "\\caption{Reduced-form sensitivity estimation at the true structural correlation: $n=300$, 500 attempts per cell. Coverage and precision use valid fits; Success uses all attempts.}",
  "\\label{tab:confounding-simulation}",
  "\\begin{tabular}{lrllrrrrr}\\toprule",
  "Error & $\\rho$ & Effect & Method & Bias & RMSE & Coverage (MCSE) & Length & Success\\\\\\midrule")
previous <- ""
for (j in seq_len(nrow(selected))) {
  d <- selected[j, ]
  if (nzchar(previous) && previous != d$Scenario) lines <- c(lines, "\\midrule")
  previous <- d$Scenario
  lines <- c(lines, sprintf("%s & %.1f & %s & %s & %.3f & %.3f & %.3f (%.3f) & %.3f & %.3f\\\\",
    ifelse(d$Scenario == "gaussian", "G", "AM"), d$TrueRho, d$Effect,
    ifelse(d$Method == "OLS", "OLS-RF", "DS-RF"), d$Bias, d$RMSE,
    d$Coverage95, d$CoverageMCSE, d$AvgLength, d$SuccessRate))
}
writeLines(c(lines, "\\bottomrule\\end{tabular}",
  "\\par\\smallskip\\begin{minipage}{0.98\\textwidth}\\footnotesize G: Gaussian; AM: asymmetric mixture. OLS-RF and DS-RF use the same reduced-form models with OLS and marginal density-score regression, respectively. They are not the conditional-outcome estimators in the original simulation. All five effects, zero-rho analyses, threshold inference and paired comparisons are retained in the project repository.\\end{minipage}",
  "\\end{table}"), "results/confounding_simulation_table.tex")

thresholds <- read_jkss_csv("results/confounding_application_thresholds.csv")
curves <- read_jkss_csv("results/confounding_application_curves.csv")
contours <- read_jkss_csv("results/confounding_application_contours.csv")
lines <- c("\\begin{table}[htbp]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
  "\\caption{JOBS II reduced-form sensitivity analysis. Nominal intervals are pointwise working-model intervals. The zero-effect correlation is conditional on a nonzero mediator treatment effect.}",
  "\\label{tab:confounding-thresholds}",
  "\\begin{tabular}{llrrr}\\toprule",
  "Effect & Method & Effect at $\\rho=0$ (95\\% CI) & $\\widehat\\rho^\\dagger$ (95\\% CI) & $(\\widehat\\rho^\\dagger)^2$\\\\\\midrule")
for (j in seq_len(nrow(thresholds))) {
  d <- thresholds[j, ]
  e <- curves[curves$Method == d$Method & curves$Effect == d$Effect & abs(curves$Rho) < 1e-10, ]
  lines <- c(lines, sprintf("%s & %s & %.3f [%.3f, %.3f] & %.3f [%.3f, %.3f] & %.3f\\\\",
    d$Effect, ifelse(d$Method == "OLS", "OLS-RF", "DS-RF"), e$Estimate, e$Lower, e$Upper,
    d$RhoZero, d$Lower, d$Upper, d$ResidualR2Product))
}
writeLines(c(lines, "\\bottomrule\\end{tabular}\\end{table}"), "results/confounding_application_table.tex")

pdf("figures/confounding_rho_curves.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), mgp = c(2.4, 0.7, 0))
for (effect in c("PNIE", "TNIE")) {
  d <- curves[curves$Effect == effect, ]
  plot(NA, xlim = c(-0.8, 0.8), ylim = range(d$Lower, d$Upper),
       xlab = expression(rho), ylab = "Natural indirect effect", main = effect)
  for (method in c("OLS", "Semiparametric")) {
    z <- d[d$Method == method, ]; z <- z[order(z$Rho), ]
    color <- ifelse(method == "OLS", "#777777", "#19677B")
    polygon(c(z$Rho, rev(z$Rho)), c(z$Lower, rev(z$Upper)), border = NA,
             col = adjustcolor(color, alpha.f = 0.15))
    lines(z$Rho, z$Estimate, col = color, lwd = 1.6, lty = ifelse(method == "OLS", 2, 1))
  }
  abline(h = 0, v = 0, lty = 3, col = "gray45")
  legend("bottomleft", c("OLS-RF", "DS-RF"), col = c("#777777", "#19677B"),
         lty = c(2, 1), bty = "n", cex = 0.85)
}
dev.off()

pdf("figures/confounding_r2_contours.pdf", width = 7, height = 6.4)
par(mfrow = c(2, 2), mar = c(3.7, 3.7, 2.7, 0.8), oma = c(2.2, 0, 0, 0), mgp = c(2.3, 0.7, 0))
axis <- sort(unique(contours$R2M))
levels <- sort(unique(c(pretty(range(contours$Estimate), n = 7), 0)))
for (method in c("OLS", "Semiparametric")) for (sign in c(-1, 1)) {
  d <- contours[contours$Method == method & contours$Sign == sign, ]
  d <- d[order(d$R2Y, d$R2M), ]
  z <- matrix(d$Estimate, nrow = length(axis))
  contour(axis, axis, z, levels = levels, xlab = expression(R[M]^"2*"),
          ylab = expression(R[Y]^"2*"), labcex = 0.65,
          main = paste(ifelse(method == "OLS", "OLS-RF", "DS-RF"),
                       ifelse(sign < 0, "opposite signs", "same signs")), cex.main = 0.9)
  if (min(z) < 0 && max(z) > 0) contour(axis, axis, z, levels = 0, add = TRUE, lwd = 2)
  for (limit in c("Lower", "Upper")) {
    z <- matrix(d[[limit]], nrow = length(axis))
    if (min(z) < 0 && max(z) > 0) contour(axis, axis, z, levels = 0, add = TRUE,
      drawlabels = FALSE, col = "#B44E28", lty = 2, lwd = 1.5)
  }
}
mtext("Black: PNIE estimate; thick black: estimate = 0; dashed orange: pointwise CI endpoint = 0",
      side = 1, outer = TRUE, cex = 0.7)
dev.off()
cat("Confounding tables and both sensitivity figures regenerated.\n")

source(file.path("R", "results_io.R"))
out <- Sys.getenv("JKSS_CURVE_OUT", "results")
figures <- Sys.getenv("JKSS_CURVE_FIGURES", "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
summary <- read_jkss_csv(file.path(out, "confounding_curve_summary.csv"))
paired <- read_jkss_csv(file.path(out, "confounding_curve_paired.csv"))
threshold <- read_jkss_csv(file.path(out, "confounding_curve_threshold_summary.csv"))
densities <- c("gaussian", "asymmetric_mixture")
correlations <- c(-0.3, 0, 0.3)
method_colors <- c(OLS = "#737373", Semiparametric = "#19677B")
density_colors <- c(gaussian = "#19677B", asymmetric_mixture = "#B44E28")
panel_title <- function(label, rho) bquote(.(label) * ": " * rho[true] == .(rho))
setup <- function() par(mfrow = c(2, 3), mar = c(3.4, 3.6, 2.8, 0.7),
                        oma = c(3.2, 0, 0, 0), mgp = c(2.2, 0.65, 0), cex = 0.84)
full_legend <- function(labels, colors, lty) {
  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", labels, col = colors, lty = lty, lwd = 1.6,
         horiz = TRUE, bty = "n", cex = 0.8)
}

pdf(file.path(figures, "confounding_simulation_curves.pdf"), width = 8.3, height = 6.4)
setup()
pnie <- summary[summary$Effect == "PNIE", ]
limits <- range(pnie$Q05, pnie$Q95, pnie$CurveTruth)
for (density in densities) for (r in correlations) {
  d <- pnie[pnie$Scenario == density & pnie$TrueRho == r, ]
  plot(NA, xlim = c(-0.8, 0.8), ylim = limits, xlab = expression(rho~"(assumed)"),
       ylab = "PNIE", main = panel_title(ifelse(density == "gaussian", "Gaussian", "Mixture"), r))
  for (method in c("OLS", "Semiparametric")) {
    z <- d[d$Method == method, ]; z <- z[order(z$AssumedRho), ]
    polygon(c(z$AssumedRho, rev(z$AssumedRho)), c(z$Q05, rev(z$Q95)),
            border = NA, col = adjustcolor(method_colors[method], 0.16))
  }
  z <- d[d$Method == "OLS", ]; z <- z[order(z$AssumedRho), ]
  lines(z$AssumedRho, z$CurveTruth, col = "black", lty = 3, lwd = 2)
  abline(v = r, col = "gray55", lty = 3)
  abline(h = -0.32, col = "#B44E28", lty = 4)
  for (method in c("OLS", "Semiparametric")) {
    z <- d[d$Method == method, ]; z <- z[order(z$AssumedRho), ]
    lines(z$AssumedRho, z$MeanEstimate, col = method_colors[method],
          lty = ifelse(method == "OLS", 2, 1), lwd = 1.4)
  }
}
full_legend(c("OLS-RF mean", "DS-RF mean", "Population curve", "Generating effect"),
            c(method_colors, "black", "#B44E28"), c(2, 1, 3, 4))
dev.off()

pdf(file.path(figures, "confounding_simulation_curve_coverage.pdf"), width = 8.3, height = 6.4)
setup()
selected <- summary[summary$Effect %in% c("PNIE", "TNIE"), ]
limits <- c(min(0.8, floor(min(selected$Coverage95) * 20) / 20), 1)
for (effect in c("PNIE", "TNIE")) for (r in correlations) {
  d <- selected[selected$Effect == effect & selected$TrueRho == r, ]
  plot(NA, xlim = c(-0.8, 0.8), ylim = limits, xlab = expression(rho~"(assumed)"),
       ylab = "Pointwise curve coverage", main = panel_title(effect, r))
  abline(h = 0.95, v = r, col = "gray65", lty = 3)
  for (density in densities) for (method in c("OLS", "Semiparametric")) {
    z <- d[d$Scenario == density & d$Method == method, ]; z <- z[order(z$AssumedRho), ]
    lines(z$AssumedRho, z$Coverage95, col = density_colors[density],
          lty = ifelse(method == "OLS", 2, 1), lwd = 1.5)
  }
}
full_legend(c("G: OLS-RF", "G: DS-RF", "AM: OLS-RF", "AM: DS-RF"),
            rep(density_colors, each = 2), c(2, 1, 2, 1))
dev.off()

pdf(file.path(figures, "confounding_simulation_curve_ratios.pdf"), width = 8.3, height = 6.4)
setup()
selected <- paired[paired$Effect %in% c("PNIE", "TNIE"), ]
limits <- range(c(1, selected$RMSERatio))
limits <- limits + c(-0.05, 0.05)
for (effect in c("PNIE", "TNIE")) for (r in correlations) {
  d <- selected[selected$Effect == effect & selected$TrueRho == r, ]
  plot(NA, xlim = c(-0.8, 0.8), ylim = limits, xlab = expression(rho~"(assumed)"),
       ylab = "RMSE ratio (DS-RF / OLS-RF)", main = panel_title(effect, r))
  abline(h = 1, v = r, col = "gray65", lty = 3)
  for (density in densities) {
    z <- d[d$Scenario == density, ]; z <- z[order(z$AssumedRho), ]
    lines(z$AssumedRho, z$RMSERatio, col = density_colors[density],
          lty = ifelse(density == "gaussian", 2, 1), lwd = 1.7)
  }
}
full_legend(c("Gaussian", "Asymmetric mixture"), density_colors, c(2, 1))
dev.off()

threshold <- threshold[order(threshold$Scenario, threshold$TrueRho, threshold$Effect, threshold$Method), ]
lines <- c("\\begin{table}[htbp]\\centering\\small\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Zero-boundary inference in the confounding simulation. The same retained data sets are used; these are not additional independent replications. Coverage and length refer to the Fisher-transform interval for $\\rho_t^\\dagger$.}",
  "\\label{tab:curve-thresholds}",
  "\\begin{tabular}{lrllrrrr}\\toprule",
  "Error & $\\rho_{\\rm true}$ & Effect & Method & Bias & RMSE & Coverage (MCSE) & Length\\\\\\midrule")
previous <- ""
for (j in seq_len(nrow(threshold))) {
  d <- threshold[j, ]
  if (nzchar(previous) && d$Scenario != previous) lines <- c(lines, "\\midrule")
  previous <- d$Scenario
  lines <- c(lines, sprintf("%s & %.1f & %s & %s & %.3f & %.3f & %.3f (%.3f) & %.3f\\\\",
    ifelse(d$Scenario == "gaussian", "G", "AM"), d$TrueRho, d$Effect,
    ifelse(d$Method == "OLS", "OLS-RF", "DS-RF"), d$Bias, d$RMSE,
    d$Coverage95, d$CoverageMCSE, d$AvgLength))
}
writeLines(c(lines, "\\bottomrule\\end{tabular}",
  sprintf("\\par\\smallskip\\begin{minipage}{0.98\\textwidth}\\footnotesize G: Gaussian; AM: asymmetric mixture. Each cell starts with %d attempts. Success denominators and failure records are the same as for Table~\\ref{tab:confounding-simulation}. A zero boundary of the point-estimate curve is not an interval-zero boundary.\\end{minipage}", unique(threshold$Attempted)),
  "\\end{table}"), file.path(out, "confounding_curve_threshold_table.tex"))
cat("Simulation sensitivity curves, coverage, paired-ratio figures and zero-boundary table written.\n")

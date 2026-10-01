source(file.path("R", "results_io.R"))
summary <- read_jkss_csv("results/calibrated_power_summary.csv")
critical <- read_jkss_csv("results/calibrated_power_critical.csv")
paired <- read_jkss_csv("results/calibrated_power_paired.csv")
methods <- c("OLS", "Huber-FIX", "Huber-SEL", "Semiparametric")
labels <- c("OLS", "Huber-FIX", "Huber-SEL", "Proposed")
format3 <- function(x) sprintf("%.3f", x)
format_interval <- function(a, b) paste0("[", format3(a), ", ", format3(b), "]")

pdf("figures/calibrated_power_curves.pdf", width = 7.2, height = 4.4,
    family = "Helvetica", useDingbats = FALSE)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), heights = c(1, 0.20))
colors <- c("#333333", "#B65C20", "#337D75", "#17518A")
shapes <- c(1, 2, 5, 16)
for (panel in seq_len(2L)) {
  par(mar = c(4.3, 4.2, 2.1, 0.7), mgp = c(2.6, 0.8, 0), las = 1)
  plot(NA, xlim = c(-0.002, 0.108), ylim = c(0, 1.015), xaxs = "i", yaxs = "i",
       xlab = expression(abs(PNIE)), ylab = "Rejection probability", xaxt = "n",
       main = if (panel == 1L) "(a) Nominal cutoff" else "(b) Null-calibrated cutoff",
       cex.main = 0.88, cex.lab = 0.9)
  axis(1, at = c(0, 0.026, 0.052, 0.078, 0.104),
       labels = c("0", ".026", ".052", ".078", ".104"), cex.axis = 0.8)
  abline(h = 0.05, lty = 3, col = "gray60")
  for (j in seq_along(methods)) {
    x <- summary[summary$Method == methods[j], ]
    x <- x[order(x$Beta2), ]
    y <- if (panel == 1L) x$NominalRate else x$CalibratedRate
    lines(abs(x$TruePNIE), y, col = colors[j], lty = j, lwd = 1.5)
    if (panel == 2L) {
      segments(abs(x$TruePNIE), x$TwoStageLower, abs(x$TruePNIE), x$TwoStageUpper,
               col = colors[j], lwd = 0.8)
    }
    points(abs(x$TruePNIE), y, col = colors[j], pch = shapes[j], cex = 0.8)
  }
}
par(mar = rep(0, 4))
plot.new()
legend("center", legend = labels, col = colors, pch = shapes, lty = seq_along(methods),
       lwd = 1.5, horiz = TRUE, bty = "n", cex = 0.85)
dev.off()

calibration_lines <- c(
  "\\begin{table}[htbp]", "\\centering\\small",
  "\\caption{Independent null-calibration results from 5,000 attempted data sets per method. The interval for the critical value reflects 1,000 resamples of the calibration tasks.}",
  "\\label{tab:power-critical}", "\\begin{tabular}{lrrr}", "\\toprule",
  "Method & Valid & Critical value & Monte Carlo 95\\% interval \\\\ \\midrule")
for (j in seq_along(methods)) {
  x <- critical[critical$Method == methods[j], ]
  calibration_lines <- c(calibration_lines, paste0(labels[j], " & ", x$Valid, " & ",
    format3(x$Critical), " & ", format_interval(x$CriticalLower, x$CriticalUpper), " \\\\"))
}
calibration_lines <- c(calibration_lines, "\\bottomrule\\end{tabular}", "\\end{table}")
writeLines(calibration_lines, "results/calibrated_power_critical_table.tex")

full <- c("\\begin{table}[htbp]", "\\centering\\small",
  "\\caption{Held-out PNIE rejection probabilities. Each effect uses 1,000 new data sets. Nominal and calibrated rates condition on valid fits; Report/reject uses all attempts. MCSE includes both calibration and evaluation uncertainty. At $\\beta_2=0$ the rates estimate type~I error, and elsewhere they estimate power.}",
  "\\label{tab:power-calibrated}", "\\begin{tabular}{rlrrrrr}", "\\toprule",
  "$\\beta_2$ & Method & Valid & Nominal & Calibrated & MCSE & Report/reject \\\\ \\midrule")
for (b in sort(unique(summary$Beta2))) {
  if (b > 0) full <- c(full, "\\midrule")
  for (j in seq_along(methods)) {
    x <- summary[summary$Beta2 == b & summary$Method == methods[j], ]
    full <- c(full, paste0(sprintf("%.1f", b), " & ", labels[j], " & ", x$Valid, " & ",
      paste(format3(c(x$NominalRate, x$CalibratedRate, x$TwoStageMCSE,
                      x$CalibratedReportReject)), collapse = " & "), " \\\\"))
  }
}
full <- c(full, "\\bottomrule\\end{tabular}", "\\end{table}")
writeLines(full, "results/calibrated_power_table.tex")

paired_lines <- c("\\begin{table}[htbp]", "\\centering\\small",
  "\\caption{Calibrated rejection contrasts on data sets where all four methods succeed. Difference is proposed minus comparator. The Monte Carlo interval uses 1,000 paired, two-stage bootstrap resamples and describes simulation uncertainty, not uncertainty in an individual mediation analysis.}",
  "\\label{tab:power-paired}", "\\begin{tabular}{rlrrrr}", "\\toprule",
  "$\\beta_2$ & Comparator & Common valid & Difference & MCSE & 95\\% interval \\\\ \\midrule")
for (b in sort(unique(paired$Beta2))) {
  if (b > 0) paired_lines <- c(paired_lines, "\\midrule")
  for (method in methods[1:3]) {
    x <- paired[paired$Beta2 == b & paired$Comparator == method, ]
    paired_lines <- c(paired_lines, paste0(sprintf("%.1f", b), " & ", method, " & ",
      x$CommonValid, " & ", format3(x$Difference), " & ", format3(x$TwoStageMCSE),
      " & ", format_interval(x$TwoStageLower, x$TwoStageUpper), " \\\\"))
  }
}
paired_lines <- c(paired_lines, "\\bottomrule\\end{tabular}", "\\end{table}")
writeLines(paired_lines, "results/calibrated_power_paired_table.tex")

macros <- character()
ids <- c("OLS", "HuberFix", "HuberSel", "Semi")
for (j in seq_along(methods)) {
  for (b in c(0, 0.1)) {
    x <- summary[summary$Beta2 == b & summary$Method == methods[j], ]
    prefix <- paste0("Cal", ids[j], if (b == 0) "Null" else "Power")
    macros <- c(macros, paste0("\\newcommand{\\", prefix, "}{", format3(x$CalibratedRate), "}"))
  }
}
sem <- summary[summary$Method == "Semiparametric", ]
weak <- sem[sem$Beta2 == 0.1, ]
contrast <- paired[paired$Beta2 == 0.1 & paired$Comparator == "Huber-SEL", ]
macros <- c(macros,
  paste0("\\newcommand{\\CalSemiSuccessRange}{", paste(format3(range(sem$SuccessRate)), collapse = "--"), "}"),
  paste0("\\newcommand{\\CalSemiReportPower}{", format3(weak$CalibratedReportReject), "}"),
  paste0("\\newcommand{\\CalHuberSelContrast}{", format3(contrast$Difference), "}"),
  paste0("\\newcommand{\\CalHuberSelContrastCI}{", format_interval(contrast$TwoStageLower, contrast$TwoStageUpper), "}"))
writeLines(macros, "results/calibrated_power_numbers.tex")
message("Calibrated-power figure, tables and manuscript macros generated.")

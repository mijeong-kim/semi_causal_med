read_effect_table_jci <- function(file, method_label, app_label) {
  dat <- read.csv(file, check.names = FALSE)
  dat$Method <- method_label
  dat$Application <- app_label
  dat
}

make_output_dir_jci <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

plot_application_effects_jci <- function(
  outfile = "figures/jci_application_effects.pdf"
) {
  jobs_ols <- read_effect_table_jci(
    "results/jci_jobs_interaction_ols_effects.csv",
    "OLS",
    "jobs"
  )
  jobs_semi <- read_effect_table_jci(
    "results/jci_jobs_interaction_semiparametric_effects.csv",
    "Semiparametric",
    "jobs"
  )
  uis_ols <- read_effect_table_jci(
    "results/jci_uis_interaction_ols_effects.csv",
    "OLS",
    "uis"
  )
  uis_semi <- read_effect_table_jci(
    "results/jci_uis_interaction_semiparametric_effects.csv",
    "Semiparametric",
    "uis"
  )

  eff_df <- rbind(jobs_ols, jobs_semi, uis_ols, uis_semi)
  effect_levels <- c("ACME(0)", "ACME(1)", "ADE(0)", "ADE(1)", "ATE")
  y_base <- rev(seq_along(effect_levels))
  names(y_base) <- effect_levels
  offset <- c(OLS = 0.15, Semiparametric = -0.15)
  ltys <- c(OLS = 1, Semiparametric = 2)
  pchs <- c(OLS = 1, Semiparametric = 16)
  panel_titles <- c(jobs = "jobs data", uis = "uis data")

  make_output_dir_jci(dirname(outfile))
  grDevices::pdf(outfile, width = 10.2, height = 5.8)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(10, 1.8))
  graphics::par(mar = c(4.8, 7.4, 3.1, 1.2), oma = c(0, 0, 0, 0))

  for (app in c("jobs", "uis")) {
    sub_df <- eff_df[eff_df$Application == app, , drop = FALSE]
    xlim <- range(c(sub_df$Lower, sub_df$Upper))
    x_pad <- 0.06 * diff(xlim)
    if (!is.finite(x_pad) || x_pad <= 0) {
      x_pad <- 0.1
    }
    plot(
      c(xlim[1] - x_pad, xlim[2] + x_pad),
      c(0.5, length(effect_levels) + 0.5),
      type = "n",
      yaxt = "n",
      ylab = "",
      xlab = "Estimate with 95% confidence interval",
      bty = "n",
      main = panel_titles[app]
    )
    graphics::axis(2, at = y_base, labels = effect_levels, las = 1)
    graphics::abline(v = 0, lty = 2, col = "gray70")
    graphics::abline(h = y_base, col = "gray92")

    for (meth in c("OLS", "Semiparametric")) {
      one <- sub_df[sub_df$Method == meth, , drop = FALSE]
      y <- y_base[one$Effect] + offset[meth]
      graphics::segments(one$Lower, y, one$Upper, y, lty = ltys[meth], lwd = 2, col = "black")
      graphics::points(one$Estimate, y, pch = pchs[meth], col = "black", cex = 1.05)
    }
  }

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c("OLS", "Semiparametric"),
    col = "black",
    pch = c(1, 16),
    lty = c(1, 2),
    lwd = 2,
    horiz = TRUE,
    bty = "n",
    cex = 0.95,
    x.intersp = 1.4,
    seg.len = 2.6
  )
}

if (sys.nframe() == 0) {
  plot_application_effects_jci()
}

source(file.path("R", "results_io.R"))
# Compare methods on identical retained data sets without changing the simulation.
source(file.path("R", "semiparametric_mediation.R"))
truth <- true_mediation_effects()
effects <- names(truth)
summaries <- list()
for (prefix in c("main_simulation", "comparator")) {
  records <- read_jkss_csv(file.path("results", paste0(prefix, "_records.csv")))
  status <- read_jkss_csv(file.path("results", paste0(prefix, "_status.csv")))
  keys <- c("Study", "SampleSize", "Scenario", "Replication", "Effect")
  cols <- c(keys, "Estimate", "Lower", "Upper")
  paired <- merge(records[records$Method == "OLS", cols],
                  records[records$Method == "Semiparametric", cols],
                  by = keys, suffixes = c("_OLS", "_Proposed"))
  stopifnot(!anyDuplicated(paired[keys]))
  groups <- split(paired, interaction(paired$SampleSize, paired$Scenario,
                                      paired$Effect, drop = TRUE))
  summaries[[prefix]] <- do.call(rbind, lapply(groups, function(g) {
    k <- g[1L, c("Study", "SampleSize", "Scenario", "Effect")]
    attempts <- status[status$SampleSize == k$SampleSize &
                       status$Scenario == k$Scenario & status$Method == "Semiparametric", ]
    stopifnot(nrow(g) == sum(attempts$Success))
    target <- unname(truth[k$Effect])
    sq_o <- (g$Estimate_OLS - target)^2
    sq_p <- (g$Estimate_Proposed - target)^2
    ratio <- sqrt(mean(sq_p) / mean(sq_o))
    ratio_mcse <- ratio * sd((sq_p / mean(sq_p) - sq_o / mean(sq_o)) / 2) / sqrt(nrow(g))
    covers_p <- g$Lower_Proposed <= target & target <= g$Upper_Proposed
    covers_o <- g$Lower_OLS <= target & target <= g$Upper_OLS
    data.frame(k, Attempts = nrow(attempts), PairedFits = nrow(g),
      SuccessRate = nrow(g) / nrow(attempts),
      BiasOLS = mean(g$Estimate_OLS - target),
      BiasProposed = mean(g$Estimate_Proposed - target),
      RMSEratio = ratio, RMSEratioMCSE = ratio_mcse,
      LengthRatio = mean(g$Upper_Proposed - g$Lower_Proposed) /
                    mean(g$Upper_OLS - g$Lower_OLS),
      CoverageOLS = mean(covers_o), CoverageProposed = mean(covers_p),
      ReportAndCover = sum(covers_p) / nrow(attempts))
  }))
}
audit <- do.call(rbind, summaries)
scenario_order <- c("Gaussian", "Skew-normal", "Asymmetric mixture", "Symmetric bimodal")
audit <- audit[order(match(audit$Study, c("Main", "Comparator")), audit$SampleSize,
                     match(audit$Scenario, scenario_order), match(audit$Effect, effects)), ]
stopifnot(nrow(audit) == 60L, all(audit$ReportAndCover <= audit$SuccessRate),
          max(abs(audit$ReportAndCover - audit$SuccessRate * audit$CoverageProposed)) < 1e-12)
write.csv(audit, "results/paired_fit_summary.csv", row.names = FALSE)
main_ng <- audit[audit$Study == "Main" & audit$Scenario != "Gaussian", ]
f <- function(x) sprintf("%.3f", x)
range_tex <- function(x) paste(f(range(x)), collapse = "--")
writeLines(c(
  paste0("\\newcommand{\\PairedRMSERange}{", range_tex(main_ng$RMSEratio), "}"),
  paste0("\\newcommand{\\PairedLengthRange}{", range_tex(main_ng$LengthRatio), "}"),
  paste0("\\newcommand{\\PairedCoverageRange}{", range_tex(main_ng$CoverageProposed), "}"),
  paste0("\\newcommand{\\ReportCoverRange}{", range_tex(main_ng$ReportAndCover), "}")
), "results/paired_fit_numbers.tex")

selected <- audit[audit$Effect %in% c("PNIE", "TE"), ]
rows <- vapply(seq_len(nrow(selected)), function(i) {
  g <- selected[i, ]
  label <- c("Gaussian" = "G", "Skew-normal" = "SN", "Asymmetric mixture" = "AM",
             "Symmetric bimodal" = "BM")[[g$Scenario]]
  paste(paste0(g$Study, " & ", g$SampleSize), label, g$Effect, g$PairedFits,
        f(g$RMSEratio), f(g$LengthRatio), f(g$CoverageOLS), f(g$CoverageProposed),
        f(g$ReportAndCover), sep = " & ")
}, character(1))
rows <- paste0(rows, " \\\\")
breaks <- which(head(selected$Scenario, -1) != tail(selected$Scenario, -1) |
                head(selected$SampleSize, -1) != tail(selected$SampleSize, -1))
rows[breaks] <- paste0(rows[breaks], "\n\\midrule")
writeLines(c(
  "\\begin{table}[htbp]\\centering\\footnotesize\\setlength{\\tabcolsep}{4pt}",
  "\\caption{Paired-fit audit using identical data sets for OLS and the proposal. Main and comparator designs each have 1,000 attempts per cell. $N$ is the number of paired valid fits. Ratios use proposed over OLS; $C_O$ and $C_P$ are conditional coverages. $R\\&C$ is the fraction of all attempts yielding a proposed interval that covers the truth, not conventional conditional coverage. G, SN, AM and BM denote the four error distributions. All five effects and ratio Monte Carlo standard errors are in the accompanying CSV file.}",
  "\\label{tab:paired-audit}",
  "\\begin{tabular}{lrllrrrrrr}\\toprule",
  "Study & $n$ & Error & Effect & $N$ & RMSE ratio & Length ratio & $C_O$ & $C_P$ & $R\\&C$ \\\\",
  "\\midrule", rows, "\\bottomrule\\end{tabular}\\end{table}"
), "results/paired_fit_table.tex")
cat("Paired-fit audit: 60 effect/design rows; no new data generated.\n")
cat("Main non-Gaussian RMSE ratios:", range_tex(main_ng$RMSEratio), "\n")
cat("Main non-Gaussian interval-length ratios:", range_tex(main_ng$LengthRatio), "\n")
cat("Report-and-cover fractions:", range_tex(main_ng$ReportAndCover), "\n")

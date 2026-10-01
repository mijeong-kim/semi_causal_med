source(file.path("R", "results_io.R"))
source(file.path("R", "compress_results.R"))

local({
  directory <- tempfile("jkss-results-io-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  path <- file.path(directory, "records.csv")
  expected <- data.frame(Estimate = rep(c(pi, NA_real_, -1e-12), 1000L),
                         Status = rep(c("valid, retained", "failed", "valid"), 1000L))
  utils::write.csv(expected, path, row.names = FALSE)
  reference <- utils::read.csv(path)
  original_hash <- unname(tools::md5sum(path))
  report <- compact_results(directory, min_bytes = 1L)
  stopifnot(nrow(report) == 1L, report$MD5 == original_hash,
            !file.exists(path), file.exists(paste0(path, ".gz")),
            jkss_result_exists(path), identical(read_jkss_csv(path), reference),
            identical(read_jkss_csv(paste0(path, ".gz")), reference),
            nrow(read_jkss_csv(path, nrows = 2L)) == 2L,
            is.null(compact_results(directory, min_bytes = 1L)))

  # A rerun's plain CSV must take precedence over an older compressed snapshot.
  reference$Estimate[1L] <- 42
  utils::write.csv(reference, path, row.names = FALSE)
  stopifnot(identical(read_jkss_csv(path), reference))
  compact_results(directory, min_bytes = 1L)
  stopifnot(!file.exists(path), identical(read_jkss_csv(path), reference))
  missing <- file.path(directory, "missing.csv")
  stopifnot(!jkss_result_exists(missing),
            inherits(try(read_jkss_csv(missing), silent = TRUE), "try-error"))
})
cat("Lossless compression, direct CSV/gzip reading and rerun precedence passed.\n")

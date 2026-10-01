copy_result_bytes <- function(source, destination) {
  repeat {
    chunk <- readBin(source, what = "raw", n = 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, destination)
  }
}

gzip_result <- function(source, destination) {
  input <- file(source, open = "rb")
  on.exit(base::close(input))
  output <- gzfile(destination, open = "wb", compression = 9L)
  on.exit(base::close(output), add = TRUE)
  copy_result_bytes(input, output)
}

gunzip_result <- function(source, destination) {
  input <- gzfile(source, open = "rb")
  on.exit(base::close(input))
  output <- file(destination, open = "wb")
  on.exit(base::close(output), add = TRUE)
  copy_result_bytes(input, output)
}

compact_result_csv <- function(path) {
  destination <- paste0(path, ".gz")
  temporary <- tempfile("jkss-compress-", tmpdir = dirname(path))
  restored <- tempfile("jkss-restore-")
  on.exit(unlink(c(temporary, restored)))
  original_hash <- unname(tools::md5sum(path))
  original_bytes <- file.info(path)$size
  gzip_result(path, temporary)
  gunzip_result(temporary, restored)
  stopifnot(identical(original_hash, unname(tools::md5sum(restored))))
  if (file.info(temporary)$size >= original_bytes) return(NULL)
  if (!file.copy(temporary, destination, overwrite = TRUE)) {
    stop("Could not write compressed result: ", destination)
  }
  # Verify the published copy and an unchanged source before removing the CSV.
  gunzip_result(destination, restored)
  stopifnot(identical(original_hash, unname(tools::md5sum(restored))),
            identical(original_hash, unname(tools::md5sum(path))))
  if (!file.remove(path)) stop("Could not remove verified redundant CSV: ", path)
  data.frame(File = basename(destination), OriginalBytes = original_bytes,
             CompressedBytes = file.info(destination)$size, MD5 = original_hash)
}

compact_results <- function(directory = "results", min_bytes = 512 * 1024) {
  if (!dir.exists(directory)) stop("Missing results directory: ", directory)
  paths <- list.files(directory, pattern = "[.]csv$", full.names = TRUE)
  paths <- paths[file.info(paths)$size >= min_bytes]
  reports <- lapply(paths, compact_result_csv)
  report <- do.call(rbind, reports)
  if (is.null(report)) {
    cat("No large uncompressed CSV files to compact.\n")
  } else {
    print(report, row.names = FALSE)
    cat(sprintf("Lossless CSV compression: %.2f MiB -> %.2f MiB.\n",
                sum(report$OriginalBytes) / 1024^2,
                sum(report$CompressedBytes) / 1024^2))
  }
  invisible(report)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) > 1L) stop("Usage: Rscript R/compress_results.R [directory]")
  compact_results(if (length(arguments)) arguments[[1L]] else "results")
}

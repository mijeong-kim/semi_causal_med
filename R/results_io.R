# Prefer newly generated CSVs; otherwise read the lossless upload copy directly.
jkss_result_exists <- function(path) {
  file.exists(path) | file.exists(paste0(path, ".gz"))
}

read_jkss_csv <- function(file, ...) {
  stopifnot(is.character(file), length(file) == 1L, !is.na(file))
  path <- if (file.exists(file)) file else paste0(file, ".gz")
  if (!file.exists(path)) {
    stop("Missing result file: ", file, " (or its .gz copy).", call. = FALSE)
  }
  if (endsWith(path, ".gz")) {
    connection <- gzfile(path, open = "rt")
    on.exit(base::close(connection))
    utils::read.csv(connection, ...)
  } else {
    utils::read.csv(path, ...)
  }
}

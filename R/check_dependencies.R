packages <- c("Matrix", "rootSolve", "sn", "mediation", "medflex", "quantreg")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
versions <- vapply(packages, function(package) {
  if (available[[package]]) as.character(utils::packageVersion(package)) else "MISSING"
}, character(1))
print(data.frame(Package = packages, Version = unname(versions)), row.names = FALSE)

if (getRversion() < "4.3.0") {
  stop("R 4.3 or later is required.", call. = FALSE)
}
if (any(!available)) {
  missing <- packages[!available]
  stop(
    "Missing or unloadable packages: ", paste(missing, collapse = ", "),
    ". Install them before running the workflow: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))",
    call. = FALSE
  )
}
cat("All R dependencies are available.\n")
cat("PDF rebuilds also require latexmk or pdflatex; ZIP creation requires zip.\n")

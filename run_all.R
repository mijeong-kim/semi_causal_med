project_files <- c(
  file.path("R", "semiparametric_mediation.R"),
  file.path("R", "run_simulations.R"),
  file.path("R", "run_application.R"),
  file.path("R", "validate_outputs.R"),
  file.path("R", "make_manuscript_assets.R"),
  "manuscript.tex",
  "supplement.tex"
)
missing_files <- project_files[!file.exists(project_files)]
if (length(missing_files)) {
  stop(
    "Run this script from the Latex/CS directory. Missing: ",
    paste(missing_files, collapse = ", ")
  )
}

run_r_script <- function(path) {
  message("\n==> ", path)
  status <- system2(file.path(R.home("bin"), "Rscript"), path)
  if (!identical(status, 0L)) {
    stop(path, " failed with status ", status)
  }
}

if (Sys.getenv("CS_SKIP_SIMULATIONS", unset = "0") != "1") {
  run_r_script(file.path("R", "run_simulations.R"))
} else {
  message("Skipping simulations because CS_SKIP_SIMULATIONS=1.")
}
run_r_script(file.path("R", "run_application.R"))
run_r_script(file.path("R", "validate_outputs.R"))
run_r_script(file.path("R", "make_manuscript_assets.R"))

dir.create(file.path("tmp", "pdfs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "pdf"), recursive = TRUE, showWarnings = FALSE)

compile_latex <- function(tex_file, output_name) {
  latexmk <- Sys.which("latexmk")
  pdflatex <- Sys.which("pdflatex")
  if (nzchar(latexmk)) {
    status <- system2(
      latexmk,
      c(
        "-pdf",
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-outdir=tmp/pdfs",
        tex_file
      )
    )
  } else if (nzchar(pdflatex)) {
    args <- c(
      "-interaction=nonstopmode",
      "-halt-on-error",
      "-output-directory=tmp/pdfs",
      tex_file
    )
    status <- 0L
    for (iteration in seq_len(3L)) {
      status <- system2(pdflatex, args)
      if (!identical(status, 0L)) break
    }
  } else {
    warning("No latexmk or pdflatex executable was found; PDFs were not compiled.")
    return(invisible(FALSE))
  }

  if (!identical(status, 0L)) {
    stop("LaTeX compilation failed for ", tex_file)
  }
  built_pdf <- file.path(
    "tmp", "pdfs", paste0(tools::file_path_sans_ext(basename(tex_file)), ".pdf")
  )
  copied <- file.copy(
    built_pdf,
    file.path("output", "pdf", output_name),
    overwrite = TRUE
  )
  if (!copied) stop("Could not copy compiled PDF: ", built_pdf)
  invisible(TRUE)
}

compile_latex("manuscript.tex", "CS_manuscript.pdf")
compile_latex("supplement.tex", "CS_supplement.pdf")

archive_files <- c(
  "README.md",
  "Makefile",
  "run_all.R",
  "manuscript.tex",
  "supplement.tex",
  list.files("R", full.names = TRUE),
  list.files("results", full.names = TRUE),
  list.files("figures", full.names = TRUE)
)
archive_files <- sort(archive_files[file.info(archive_files)$isdir %in% FALSE])
archive_path <- file.path(
  "output", "CS_Online_Resource_2_reproducibility.zip"
)
if (file.exists(archive_path) && !file.remove(archive_path)) {
  stop("Could not replace the existing reproducibility archive.")
}
utils::zip(archive_path, archive_files, flags = "-r9X")
if (!file.exists(archive_path) || file.info(archive_path)$size == 0) {
  stop("The reproducibility archive was not created.")
}

message("\nComputational Statistics reproduction workflow and Online Resource 2 archive completed.")

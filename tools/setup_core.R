# Run with R >= 4.3 and GRSEM_PROJECT_ROOT set; installs only Core/dev dependencies.
root <- Sys.getenv("GRSEM_PROJECT_ROOT", "")
if (!nzchar(root)) stop("set GRSEM_PROJECT_ROOT")
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
stopifnot(unname(read.dcf(file.path(root, "DESCRIPTION"))[1, "Package"]) == "grsem")
lib <- file.path(root, ".Rlib")
downloads <- file.path(root, "artifacts", "downloads")
dir.create(lib, showWarnings = FALSE)
dir.create(downloads, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
install.packages(c("lavaan", "testthat", "numDeriv", "pkgload", "rcmdcheck",
                   "knitr", "rmarkdown"),
                 lib = lib, repos = "https://cloud.r-project.org",
                 dependencies = NA, destdir = downloads)
if (packageVersion("lavaan") < "0.7.2") {
  stop("lavaan >= 0.7-2 with lav_export_estimation() is required; do not substitute an older backend")
}

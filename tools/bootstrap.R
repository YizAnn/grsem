# For repository scripts only; never installs dependencies or reads user data.
root <- Sys.getenv("GRSEM_PROJECT_ROOT", "")
if (!nzchar(root)) stop("set GRSEM_PROJECT_ROOT to this candidate repository")
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
options(grsem.project_root = root)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(grsem))
loaded <- normalizePath(find.package("grsem"), winslash = "/", mustWork = TRUE)
if (!startsWith(tolower(loaded), paste0(tolower(root), "/.rlib/"))) {
  stop("install the candidate package in its .Rlib before running repository scripts")
}
output_dir <- function(relative) {
  path <- grsem:::.grsem_output_path(file.path("artifacts", relative, ".anchor"),
                                    create_parent = TRUE)
  dirname(path)
}

# Validate each output file too, so an existing symlink cannot escape the root.
write.csv <- function(x, file, ...) {
  utils::write.csv(x, grsem:::.grsem_output_path(file, TRUE), ...)
}
writeLines <- function(text, con, ...) {
  base::writeLines(text, grsem:::.grsem_output_path(con, TRUE), ...)
}
saveRDS <- function(object, file, ...) {
  base::saveRDS(object, grsem:::.grsem_output_path(file, TRUE), ...)
}

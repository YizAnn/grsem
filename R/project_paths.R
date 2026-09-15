.grsem_project_root <- function() {
  root <- getOption("grsem.project_root", Sys.getenv("GRSEM_PROJECT_ROOT", ""))
  if (!nzchar(root)) {
    root <- getwd()
    repeat {
      desc <- file.path(root, "DESCRIPTION")
      if (file.exists(desc) &&
          identical(unname(read.dcf(desc)[1, "Package"]), "grsem")) break
      parent <- dirname(root)
      if (identical(parent, root)) {
        stop("set GRSEM_PROJECT_ROOT to the candidate repository for file outputs")
      }
      root <- parent
    }
  }
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  desc <- file.path(root, "DESCRIPTION")
  if (!file.exists(desc) ||
      !identical(unname(read.dcf(desc)[1, "Package"]), "grsem")) {
    stop("project root must be a grsem repository")
  }
  root
}

.grsem_output_path <- function(path, create_parent = FALSE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) stop("invalid output path")
  root <- .grsem_project_root()
  path <- gsub("\\\\", "/", path)
  if (any(strsplit(path, "/", fixed = TRUE)[[1L]] == "..")) {
    stop("output path must not contain parent traversal")
  }
  if (!grepl("^(/|[A-Za-z]:)", path)) path <- file.path(root, path)
  # Resolve existing ancestors before creating anything (including symlinks).
  ancestor <- path
  suffix <- character()
  while (!file.exists(ancestor) && !dir.exists(ancestor)) {
    suffix <- c(basename(ancestor), suffix)
    parent <- dirname(ancestor)
    if (identical(parent, ancestor)) stop("cannot resolve output path")
    ancestor <- parent
  }
  resolved <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
  if (length(suffix)) resolved <- paste(c(resolved, suffix), collapse = "/")
  a <- if (.Platform$OS.type == "windows") tolower(resolved) else resolved
  b <- if (.Platform$OS.type == "windows") tolower(root) else root
  if (!(identical(a, b) || startsWith(a, paste0(b, "/")))) {
    stop("output path must stay inside the candidate project")
  }
  if (create_parent) {
    dir.create(dirname(resolved), recursive = TRUE, showWarnings = FALSE)
    # Re-check after creation before the caller writes.
    return(.grsem_output_path(resolved, create_parent = FALSE))
  }
  resolved
}

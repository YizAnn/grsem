# Stage only R-package inputs so an in-repository tempdir cannot copy itself.
source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
stage_parent <- output_dir(paste0("build-source-", format(Sys.time(), "%Y%m%d%H%M%S"),
                                  "-", Sys.getpid()))
stage <- file.path(stage_parent, "grsem")
if (dir.exists(stage)) stop("staging directory already exists")
dir.create(stage)
inputs <- c("DESCRIPTION", "NAMESPACE", "LICENSE", "NEWS.md", "README.md",
            ".Rbuildignore", "R", "man", "tests", "vignettes", "inst")
for (item in inputs) {
  if (!file.copy(file.path(root, item), stage, recursive = TRUE,
                 copy.mode = TRUE, copy.date = TRUE)) stop("staging failed: ", item)
}
writeLines(stage, file.path(root, "artifacts", "build-source-path.txt"))
cat("Staged source:", stage, "\n")

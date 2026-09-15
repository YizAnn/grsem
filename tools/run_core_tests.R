source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
library(testthat)
result <- testthat::test_local(root, reporter = "summary", stop_on_failure = FALSE)
table <- as.data.frame(result)
write.csv(table[, setdiff(names(table), "result"), drop = FALSE],
          grsem:::.grsem_output_path("artifacts/core-tests.csv", TRUE), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           grsem:::.grsem_output_path("artifacts/sessionInfo.txt", TRUE))
if (any(table$failed > 0 | table$error)) quit(status = 1L)

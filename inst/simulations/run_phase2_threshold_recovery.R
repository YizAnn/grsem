source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
## Predeclared Phase 2 categorical recovery supplement.
## It is separate from the core job because that process loaded its source
## before threshold RMSE was added to the result schema.
suppressPackageStartupMessages(library(lavaan))

B <- 100L
n <- 500L
seed0 <- 2026071902L
load_truth <- c(0.80, 0.70, 0.60, 0.75, 0.65, 0.55)
threshold_truth <- qnorm(c(0.90, 0.93, 0.95, 0.90, 0.94, 0.96))
items <- paste0("T", seq_along(load_truth))
model <- paste("F =~", paste(items, collapse = " + "))

one_rep <- function(b) {
  set.seed(seed0 + b)
  eta <- rnorm(n)
  latent <- vapply(seq_along(items), function(j) {
    load_truth[j] * eta + sqrt(1 - load_truth[j]^2) * rnorm(n)
  }, numeric(n))
  dat <- as.data.frame(sweep(latent, 2, threshold_truth, ">"))
  names(dat) <- items
  start <- proc.time()[[3L]]
  fit <- try(cfa(model, dat, ordered = items, estimator = "DWLS",
                 parameterization = "delta", std.lv = TRUE,
                 se = "none", test = "none"), silent = TRUE)
  if (inherits(fit, "try-error") || !lavInspect(fit, "converged")) {
    return(data.frame(replicate = b, status = "failed",
      loading_rmse = NA_real_, threshold_rmse = NA_real_,
      runtime_seconds = proc.time()[[3L]] - start,
      message = if (inherits(fit, "try-error")) as.character(fit) else
        "lavaan did not converge"))
  }
  pt <- parameterEstimates(fit)
  loads <- setNames(pt$est[pt$op == "=~"], pt$rhs[pt$op == "=~"])[items]
  thres <- setNames(pt$est[pt$op == "|"], pt$lhs[pt$op == "|"])[items]
  data.frame(replicate = b, status = "converged",
    loading_rmse = sqrt(mean((loads - load_truth)^2)),
    threshold_rmse = sqrt(mean((thres - threshold_truth)^2)),
    runtime_seconds = proc.time()[[3L]] - start, message = "")
}

raw <- do.call(rbind, lapply(seq_len(B), one_rep))
ok <- raw$status == "converged"
summary <- data.frame(
  replications = B, failure_rate = mean(!ok),
  mean_loading_rmse = mean(raw$loading_rmse[ok], na.rm = TRUE),
  mean_threshold_rmse = mean(raw$threshold_rmse[ok], na.rm = TRUE),
  mean_runtime_seconds = mean(raw$runtime_seconds[ok], na.rm = TRUE),
  seed_base = seed0
)
out <- output_dir("simulations/results_phase2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(raw, file.path(out, "phase2_threshold_recovery_raw.csv"),
          row.names = FALSE)
write.csv(summary, file.path(out, "phase2_threshold_recovery_summary.csv"),
          row.names = FALSE)
print(summary)

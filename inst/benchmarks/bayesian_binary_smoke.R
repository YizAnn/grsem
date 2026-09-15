source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
stop("Experimental Bayesian benchmark is disabled in the Core build; no sampling or compilation is run.")
suppressPackageStartupMessages({
  library(grsem)
  library(cmdstanr)
})
# CmdStan is intentionally not configured by the Core workflow.
set.seed(20260721)
N <- 200L; P <- 4L; K <- 1L; G <- 2L
eta <- rnorm(N)
loading_true <- c(1, 0.8, 0.6, 0)
threshold_true <- c(0.2, 0.7, 1.1, 1.4)
prob <- sapply(seq_len(P), function(p)
  pnorm(loading_true[p] * eta - threshold_true[p]))
y <- matrix(rbinom(N * P, 1, as.vector(prob)), N, P)
loading_allowed <- matrix(1L, P, K)
loading_allowed[1, 1] <- 0L
loading_fixed <- matrix(0, P, K)
loading_fixed[1, 1] <- 1
loading_group <- matrix(c(0L, 1L, 1L, 2L), P, K)
stan_data <- list(
  N = N, P = P, K = K, y = y,
  loading_allowed = loading_allowed,
  loading_fixed = loading_fixed,
  loading_group = loading_group, G = G,
  global_scale = 0.2
)
spec <- grsem_model(
  "f =~ y1 + y2 + y3 + y4",
  data.frame(y1 = y[, 1], y2 = y[, 2], y3 = y[, 3], y4 = y[, 4]),
  ordered = paste0("y", 1:4), metadata = list(stan_data = stan_data)
)
fit <- fit_bayes(
  spec, priors = grsem_prior("group_local", global_scale = 0.2),
  stan_model = "binary_probit", chains = 2, parallel_chains = 2,
  seed = 20260721, iter_warmup = 200, iter_sampling = 200,
  refresh = 100, adapt_delta = 0.90
)
out <- output_dir("benchmarks")
variables <- c("loading[1,1]", "loading[2,1]", "loading[3,1]",
               "loading[4,1]", paste0("threshold[", 1:4, "]"))
summary <- fit$fit$summary(variables = variables)
summary$truth <- c(loading_true, threshold_true)
write.csv(summary, file.path(out, "bayesian_binary_smoke_summary.csv"),
          row.names = FALSE)
write.csv(fit$diagnostics,
          file.path(out, "bayesian_binary_smoke_diagnostics.csv"),
          row.names = FALSE)
saveRDS(fit, file.path(out, "bayesian_binary_smoke_fit.rds"))
print(summary)
print(fit$diagnostics)

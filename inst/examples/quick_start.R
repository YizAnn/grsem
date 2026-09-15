# Completely synthetic Core workflow. Run in a fresh, installed-package session.
library(grsem)
set.seed(20260915)
n <- 400L
x <- rnorm(n)
z <- rnorm(n)
m <- 0.8 * x + rnorm(n, sd = 0.7)
y <- 0.7 * m + 0.25 * x + rnorm(n, sd = 0.7)
dat <- data.frame(x = x, z = z, m = m, y = y)
spec <- grsem_model("m ~ x + z\ny ~ m + x + z", dat)
groups <- c("m ~ x" = "m_inputs", "m ~ z" = "m_inputs",
            "y ~ m" = "y_inputs", "y ~ x" = "y_inputs", "y ~ z" = "y_inputs")
fit <- fit_penalized(
  spec, groups = groups, penalize = "regressions",
  rho = c(0.2, 0), alpha = 0.5,
  control = grsem_control(algorithm = "pg", max_iter = 3000,
                          tolerance = 1e-5, n_starts = 1)
)
# Predeclared demonstration setting, NOT a CV-selected optimum.
solution <- fit$path$solution[fit$path$rho == 0.2]
stopifnot(length(solution) == 1L, fit$path$status[match(solution, fit$path$solution)] != "failed")
edges <- selected_edges(fit, solution)
dag <- validate_dag(edges, nodes = names(dat))
print(fit$path[, c("solution", "rho", "status", "prox_gradient_max")])
print(edges)
print(dag)
stopifnot(dag$is_dag)
refit <- post_selection_refit(fit, solution, data = dat)
stopifnot(isTRUE(lavaan::lavInspect(refit, "converged")))
cat("Synthetic Quick Start completed; graph legality is not causal identification.\n")

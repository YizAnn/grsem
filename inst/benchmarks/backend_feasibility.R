source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
# Minimal route-A/route-B feasibility benchmark.
# Run from the package root with the project R 4.4.3 wrapper.
suppressPackageStartupMessages({
  library(lavaan)
  library(numDeriv)
})

set.seed(20260718)
n <- 500L
eta <- rnorm(n)
dat <- data.frame(
  x1 = 0.9 * eta + rnorm(n, sd = 0.5),
  x2 = 0.8 * eta + rnorm(n, sd = 0.6),
  x3 = 0.7 * eta + rnorm(n, sd = 0.7),
  x4 = 0.5 * eta + rnorm(n, sd = 0.8)
)
model <- "f =~ x1 + x2 + x3 + x4"

t_lavaan <- system.time({
  fit_lav <- cfa(model, data = dat, std.lv = TRUE, se = "none", test = "none")
})[["elapsed"]]
exported <- lav_export_estimation(fit_lav)
x <- exported$starting_values + seq_along(exported$starting_values) * 1e-5
t_analytic <- system.time(g_analytic <- exported$gradient_function(x, fit_lav))[["elapsed"]]
t_numeric <- system.time(g_numeric <- numDeriv::grad(
  function(z) exported$objective_function(z, fit_lav), x
))[["elapsed"]]

rows <- data.frame(
  route = "A_lavaan_export",
  package_version = as.character(packageVersion("lavaan")),
  unpenalized_fit_seconds = t_lavaan,
  analytic_gradient_seconds = t_analytic,
  numeric_gradient_seconds = t_numeric,
  max_gradient_error = max(abs(g_analytic - g_numeric)),
  equality_constraint_handling = "lavaan K mapping; simple equalities tested",
  exact_group_prox_in_native_backend = FALSE,
  role = "selected main backend",
  stringsAsFactors = FALSE
)

if (requireNamespace("OpenMx", quietly = TRUE)) {
  suppressPackageStartupMessages(library(OpenMx))
  manifests <- names(dat)
  mx_model <- mxModel(
    "one_factor", type = "RAM", manifestVars = manifests, latentVars = "f",
    mxPath(from = "f", to = manifests, arrows = 1, free = TRUE,
           values = c(0.9, 0.8, 0.7, 0.5), labels = paste0("l", 1:4)),
    mxPath(from = manifests, arrows = 2, free = TRUE, values = 0.5),
    mxPath(from = "f", arrows = 2, free = FALSE, values = 1),
    mxPath(from = "one", to = manifests, arrows = 1,
           free = TRUE, values = 0),
    mxData(dat, type = "raw")
  )
  t_mx <- system.time(fit_mx <- mxRun(mx_model, silent = TRUE))[["elapsed"]]
  rows <- rbind(rows, data.frame(
    route = "B_OpenMx_RAM",
    package_version = as.character(packageVersion("OpenMx")),
    unpenalized_fit_seconds = t_mx,
    analytic_gradient_seconds = NA_real_,
    numeric_gradient_seconds = NA_real_,
    max_gradient_error = max(abs(fit_mx$output$gradient), na.rm = TRUE),
    equality_constraint_handling = "native constraints strong; separate RAM mapping required",
    exact_group_prox_in_native_backend = FALSE,
    role = "reference/future FIML backend",
    stringsAsFactors = FALSE
  ))
}

write.csv(rows, file.path(output_dir("benchmarks"), "backend_feasibility_results.csv"), row.names = FALSE)
print(rows)

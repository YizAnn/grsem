.grsem_cmdstan_available <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  !inherits(version, "try-error") && length(version) == 1L
}

#' Fit an experimental Bayesian SEM backend
#'
#' The Phase 2 interface deliberately delegates HMC/NUTS to CmdStan. A model
#' specification must provide `metadata$stan_data`; the research prototype does
#' not silently translate unsupported lavaan constraints.
#'
#' @param spec A [grsem_model()] specification.
#' @param priors A [grsem_prior()] specification.
#' @param backend Currently only `"cmdstanr"`.
#' @param stan_model Continuous, binary-probit, or hierarchical-version model.
#' @param chains Number of chains.
#' @param parallel_chains Parallel chains.
#' @param seed Random seed.
#' @param ... Sampling arguments passed to CmdStan.
#' @return A `grsem_bayes_fit` object.
#' @export
fit_bayes <- function(spec, priors = grsem_prior(), backend = "cmdstanr",
                      stan_model = c("continuous", "binary_probit",
                                     "hierarchical_versions"),
                      chains = 4L, parallel_chains = chains,
                      seed = 1L, ...) {
  if (!inherits(spec, "grsem_model_spec")) {
    stop("spec must be created by grsem_model()", call. = FALSE)
  }
  stan_model <- match.arg(stan_model)
  if (!identical(backend, "cmdstanr")) {
    stop("Phase 2 Bayesian estimation delegates only to cmdstanr", call. = FALSE)
  }
  if (!.grsem_cmdstan_available()) {
    stop("cmdstanr is installed but CmdStan is not configured; Bayesian estimation is not implemented in this environment",
         call. = FALSE)
  }
  stan_data <- spec$metadata$stan_data
  if (is.null(stan_data)) {
    stop("spec$metadata$stan_data is required; automatic general lavaan-to-Stan translation is not implemented",
         call. = FALSE)
  }
  file <- system.file("stan", paste0("sem_", stan_model, ".stan"),
                      package = "grsem")
  if (!nzchar(file)) {
    file <- file.path("inst", "stan", paste0("sem_", stan_model, ".stan"))
  }
  exe <- sub("[.]stan$", ".exe", file)
  model <- if (file.exists(exe)) {
    cmdstanr::cmdstan_model(stan_file = file, exe_file = exe, compile = FALSE)
  } else {
    cmdstanr::cmdstan_model(file)
  }
  fit <- model$sample(
    data = stan_data, chains = chains, parallel_chains = parallel_chains,
    seed = seed, ...
  )
  diagnostics <- fit$diagnostic_summary()
  structure(list(
    fit = fit, spec = spec, priors = priors, stan_model = stan_model,
    diagnostics = diagnostics, status = if (all(diagnostics$num_divergent == 0))
      "converged" else "approximately_converged"
  ), class = "grsem_bayes_fit")
}

#' Backward-compatible Bayesian interface
#' @inheritParams fit_bayes
#' @export
grsem_bayes <- function(spec, priors = grsem_prior(), backend = "cmdstanr",
                        ...) {
  fit_bayes(spec, priors = priors, backend = backend, ...)
}

#' @export
summary.grsem_bayes_fit <- function(object, ...) {
  list(
    posterior = object$fit$summary(),
    diagnostics = object$diagnostics,
    status = object$status,
    interpretation = "Posterior summaries are shrinkage estimates, not exact-zero selections."
  )
}

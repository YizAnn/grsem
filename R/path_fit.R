#' Bootstrap a relative regularization path
#'
#' Each resample recomputes `lambda_max` and maps the requested relative path
#' to its own absolute lambda values. The function reports failures rather than
#' dropping them.
#'
#' @param spec A [grsem_model()] specification.
#' @param groups Parameter groups passed to [fit_penalized()].
#' @param rho Relative regularization grid.
#' @param alpha Mixing grid.
#' @param B Number of bootstrap resamples.
#' @param seed Random seed.
#' @param ... Additional arguments passed to [fit_penalized()].
#' @return List with fold-level diagnostics and fitted objects.
#' @export
bootstrap_penalized <- function(spec, groups, rho, alpha = c(0, 0.5, 1),
                                B = 100L, seed = 1L, ...) {
  if (!inherits(spec, "grsem_model_spec")) {
    stop("spec must be created by grsem_model()", call. = FALSE)
  }
  set.seed(seed)
  diagnostics <- vector("list", B)
  fits <- vector("list", B)
  for (b in seq_len(B)) {
    index <- sample.int(nrow(spec$data), replace = TRUE)
    boot_spec <- spec
    boot_spec$data <- spec$data[index, , drop = FALSE]
    fit <- try(fit_penalized(
      boot_spec, groups = groups, rho = rho, alpha = alpha, ...
    ), silent = TRUE)
    if (inherits(fit, "try-error")) {
      diagnostics[[b]] <- data.frame(
        bootstrap = b, alpha = NA_real_, lambda_max = NA_real_,
        rho = NA_real_, lambda = NA_real_, nonzero_penalized = NA_integer_,
        nonzero_groups = NA_integer_, prox_gradient_max = NA_real_,
        objective = NA_real_, status = "failed",
        message = as.character(fit), stringsAsFactors = FALSE
      )
      next
    }
    rows <- fit$path[, c(
      "alpha", "lambda_max", "rho", "lambda", "nonzero_penalized",
      "nonzero_groups", "prox_gradient_max", "objective", "status", "message"
    )]
    rows$bootstrap <- b
    diagnostics[[b]] <- rows[, c("bootstrap", setdiff(names(rows), "bootstrap"))]
    fits[[b]] <- fit
  }
  structure(list(
    diagnostics = do.call(rbind, diagnostics), fits = fits,
    rho = rho, alpha = alpha, seed = seed
  ), class = "grsem_bootstrap")
}

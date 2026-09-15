#' Fit a penalized SEM from a unified specification
#'
#' @param spec A [grsem_model()] specification.
#' @param groups Parameter-group definition passed to [grsem()].
#' @param backend Continuous ML or experimental binary DWLS backend.
#' @param ... Additional arguments passed to [grsem()].
#' @return A `grsem_fit` object retaining the originating specification.
#' @export
fit_penalized <- function(spec, groups = NULL,
                          backend = if (length(spec$ordered))
                            "lavaan_dwls" else "lavaan_ml", ...) {
  if (!inherits(spec, "grsem_model_spec")) {
    stop("spec must be created by grsem_model()", call. = FALSE)
  }
  fit <- grsem(
    model = spec$model, data = spec$data, groups = groups,
    backend = backend, ordered = spec$ordered, ...
  )
  fit$model_spec <- spec
  fit
}

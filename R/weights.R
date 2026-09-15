#' Construct group weights for a grsem model
#'
#' @param model,data,groups,penalize,estimator,std.lv,meanstructure Model
#'   arguments as in [grsem()].
#' @param method One of `"sqrt_size"`, `"none"`, or `"hessian"`.
#' @param ... Additional arguments passed to lavaan.
#' @return Named positive group weights.
#' @export
grsem_group_weights <- function(model, data, groups,
                                penalize = "loadings",
                                method = c("sqrt_size", "none", "hessian"),
                                estimator = "ML", std.lv = TRUE,
                                meanstructure = FALSE, ...) {
  method <- match.arg(method)
  backend <- .grsem_prepare_backend(
    model, data, estimator, std.lv, meanstructure, ...
  )
  spec <- .grsem_build_spec(
    backend, groups, penalize,
    penalty_factor = NULL, group_weights = NULL,
    ungrouped = "unpenalized"
  )
  ids <- names(spec$group_weights)
  if (method == "sqrt_size") return(spec$group_weights)
  if (method == "none") return(stats::setNames(rep(1, length(ids)), ids))
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("numDeriv is required for Hessian group weights", call. = FALSE)
  }
  hessian <- numDeriv::hessian(backend$fn, backend$x_unpenalized)
  raw <- vapply(ids, function(id) {
    idx <- which(spec$groups == id & spec$penalty_factor > 0)
    block_diag <- pmax(diag(hessian)[idx], .Machine$double.eps)
    sqrt(sum(block_diag))
  }, numeric(1))
  if (any(!is.finite(raw)) || any(raw <= 0)) {
    stop("the unpenalized Hessian did not yield valid group weights",
         call. = FALSE)
  }
  raw / exp(mean(log(raw)))
}

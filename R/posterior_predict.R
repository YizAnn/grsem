#' Extract posterior predictive draws
#' @param object A `grsem_bayes_fit` object.
#' @param variable Generated-quantity variable name.
#' @param ... Passed to CmdStan draws extraction.
#' @export
posterior_predict <- function(object, variable = "y_rep", ...) {
  if (!inherits(object, "grsem_bayes_fit")) stop("not a grsem_bayes_fit")
  object$fit$draws(variables = variable, ...)
}

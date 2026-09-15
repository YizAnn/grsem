#' Draw prior predictive quantities from a fitted CmdStan model
#' @param object A `grsem_bayes_fit` object.
#' @param ... Passed to the CmdStan generated-quantities method.
#' @export
prior_predict <- function(object, ...) {
  if (!inherits(object, "grsem_bayes_fit")) stop("not a grsem_bayes_fit")
  stop("prior predictive execution requires a dedicated fixed_param program; interface reserved but not implemented",
       call. = FALSE)
}

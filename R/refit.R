#' Post-selection unpenalized SEM refit
#'
#' Fixes zeroed penalized parameters to zero and re-estimates the selected model
#' without a penalty. This does not make ordinary p-values unconditional.
#'
#' @param object A `grsem_fit` object.
#' @param solution Solution identifier.
#' @param data Training data used for refitting; defaults to the original data.
#' @param store Store the refit in the returned object.
#' @param ... Additional arguments passed to [lavaan::sem()].
#' @return If `store=TRUE`, an updated `grsem_fit`; otherwise a lavaan object.
#' @export
post_selection_refit <- function(object, solution, data = object$data,
                                 store = FALSE, ...) {
  stopifnot(inherits(object, "grsem_fit"))
  if (!identical(object$backend$estimator, "ML")) {
    stop("post_selection_refit currently supports the Core ML backend only", call. = FALSE)
  }
  row <- match(solution, object$path$solution)
  if (is.na(row)) stop("unknown solution", call. = FALSE)
  x <- object$coefficients[[row]]
  pt <- object$spec$full_parameter_table
  free_ids <- sort(unique(pt$free[pt$free > 0L]))
  for (i in seq_along(free_ids)) {
    opt_index <- object$spec$full_to_opt[i]
    if (opt_index == 0L) next
    if (object$spec$penalty_factor[opt_index] > 0 &&
        abs(x[opt_index]) <= object$control$zero_tolerance) {
      hit <- which(pt$free == free_ids[i])
      pt$free[hit] <- 0L
      pt$ustart[hit] <- 0
      if ("est" %in% names(pt)) pt$est[hit] <- 0
    }
  }
  refit <- lavaan::sem(
    model = pt, data = data, estimator = "ML",
    std.lv = object$backend$std.lv,
    meanstructure = object$backend$meanstructure,
    ...
  )
  attr(refit, "selection_warning") <-
    "Structure was selected on the same training data; ordinary p-values are conditional, not unconditional post-selection inference."
  if (!store) return(refit)
  object$refits[[as.character(solution)]] <- refit
  object
}

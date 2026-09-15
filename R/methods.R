#' @export
print.grsem_fit <- function(x, ...) {
  cat("Grouped regularized SEM research prototype\n")
  cat("Solutions:", nrow(x$path), "\n")
  cat("Status:\n")
  print(table(x$path$status, useNA = "ifany"))
  cat("Backend: lavaan", as.character(utils::packageVersion("lavaan")), "\n")
  invisible(x)
}

#' @export
summary.grsem_fit <- function(object, solution = NULL, ...) {
  if (is.null(solution)) return(object$path)
  row <- match(solution, object$path$solution)
  if (is.na(row)) stop("unknown solution", call. = FALSE)
  list(
    diagnostics = object$path[row, , drop = FALSE],
    parameter_table = .grsem_parameter_table(object, object$coefficients[[row]]),
    selected_groups = selected_groups(object, solution)
  )
}

#' @export
coef.grsem_fit <- function(object, solution = NULL, ...) {
  if (is.null(solution)) solution <- tail(object$path$solution, 1L)
  row <- match(solution, object$path$solution)
  if (is.na(row)) stop("unknown solution", call. = FALSE)
  if (is.function(object$backend$get_coef)) {
    object$backend$get_coef(object$coefficients[[row]])
  } else {
    object$backend$exported$get_coef(
      object$coefficients[[row]], object$backend$fit
    )
  }
}

#' Return selected parameter groups
#'
#' @param object A `grsem_fit` object.
#' @param solution Solution identifier.
#' @param tolerance Optional zero threshold.
#' @return Data frame with group IDs, norms, and selection status.
#' @export
selected_groups <- function(object, solution, tolerance = NULL) {
  stopifnot(inherits(object, "grsem_fit"))
  row <- match(solution, object$path$solution)
  if (is.na(row)) stop("unknown solution", call. = FALSE)
  if (is.null(tolerance)) tolerance <- object$control$zero_tolerance
  x <- object$coefficients[[row]]
  ids <- names(object$spec$group_weights)
  data.frame(
    group = ids,
    size = vapply(ids, function(id) sum(object$spec$groups == id, na.rm = TRUE),
                  integer(1)),
    norm = vapply(ids, function(id) sqrt(sum(x[which(object$spec$groups == id)]^2)),
                  numeric(1)),
    selected = vapply(ids, function(id) {
      any(abs(x[which(object$spec$groups == id)]) > tolerance)
    }, logical(1)),
    stringsAsFactors = FALSE
  )
}

#' @export
plot.grsem_fit <- function(x, type = c("path", "diagnostic"), ...) {
  type <- match.arg(type)
  if (type == "diagnostic") {
    graphics::plot(x$path$rho, x$path$prox_gradient_max,
                   log = "x", xlab = "rho = lambda/lambda_max",
                   ylab = "maximum proximal-gradient mapping",
                   pch = c(converged = 16, approximately_converged = 17,
                           failed = 1)[x$path$status], ...)
    return(invisible(x))
  }
  alpha_values <- unique(x$path$alpha)
  old <- graphics::par(mfrow = c(length(alpha_values), 1L))
  on.exit(graphics::par(old), add = TRUE)
  for (a in alpha_values) {
    idx <- which(x$path$alpha == a)
    mat <- do.call(rbind, x$coefficients[idx])
    pen <- which(x$spec$penalty_factor > 0)
    graphics::matplot(x$path$rho[idx], mat[, pen, drop = FALSE],
                      type = "l", lty = 1, log = "x",
                      xlab = "rho = lambda/lambda_max",
                      ylab = "penalized estimate",
                      main = paste("alpha =", a), ...)
  }
  invisible(x)
}

#' @export
predict.grsem_fit <- function(object, newdata = NULL, solution = NULL,
                              type = "lv", ...) {
  if (is.null(solution)) solution <- tail(object$path$solution, 1L)
  key <- as.character(solution)
  if (is.null(object$refits[[key]])) {
    stop("predict() is available after post_selection_refit(); penalized estimates are not treated as ordinary fitted lavaan models",
         call. = FALSE)
  }
  lavaan::lavPredict(object$refits[[key]], newdata = newdata, type = type, ...)
}

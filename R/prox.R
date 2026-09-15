#' Soft-thresholding operator
#'
#' @param x Numeric vector.
#' @param threshold Non-negative scalar or vector.
#' @return Thresholded vector.
#' @export
soft_threshold <- function(x, threshold) {
  if (any(!is.finite(threshold)) || any(threshold < 0)) {
    stop("threshold must be finite and non-negative", call. = FALSE)
  }
  sign(x) * pmax(abs(x) - threshold, 0)
}

#' Group soft-thresholding operator
#'
#' @param x Numeric group vector.
#' @param threshold Non-negative scalar.
#' @return Thresholded group vector.
#' @export
group_soft_threshold <- function(x, threshold) {
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    stop("threshold must be one finite non-negative value", call. = FALSE)
  }
  norm_x <- sqrt(sum(x * x))
  if (norm_x == 0 || norm_x <= threshold) {
    return(rep.int(0, length(x)))
  }
  (1 - threshold / norm_x) * x
}

#' Sparse-group proximal operator
#'
#' Applies elementwise soft thresholding followed by non-overlapping group
#' soft thresholding. Parameters with `penalty_factor == 0` are unchanged.
#'
#' @param x Numeric parameter vector.
#' @param groups Character group ID per coordinate; `NA` denotes no group.
#' @param step_lambda Product of step size and lambda.
#' @param alpha Mixing value in `[0, 1]`.
#' @param group_weights Named group-weight vector.
#' @param penalty_factor Per-coordinate L1 factor. The first prototype supports
#'   zero/one factors for grouped coordinates.
#' @param l1_only Logical vector identifying ungrouped L1 coordinates.
#' @return Proximal update.
#' @export
sparse_group_prox <- function(x, groups, step_lambda, alpha,
                              group_weights, penalty_factor,
                              l1_only = rep(FALSE, length(x))) {
  stopifnot(length(x) == length(groups),
            length(x) == length(penalty_factor),
            length(x) == length(l1_only))
  if (!is.finite(step_lambda) || step_lambda < 0 ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("step_lambda and alpha are outside their valid ranges", call. = FALSE)
  }
  out <- x
  penalized <- penalty_factor > 0 & (!is.na(groups) | l1_only)
  out[penalized] <- soft_threshold(
    out[penalized],
    step_lambda * alpha * penalty_factor[penalized]
  )
  group_ids <- unique(groups[!is.na(groups) & penalty_factor > 0])
  for (group_id in group_ids) {
    idx <- which(groups == group_id & penalty_factor > 0)
    weight <- unname(group_weights[group_id])
    if (length(weight) != 1L || !is.finite(weight) || weight <= 0) {
      stop("missing or invalid weight for group ", group_id, call. = FALSE)
    }
    out[idx] <- group_soft_threshold(
      out[idx], step_lambda * (1 - alpha) * weight
    )
  }
  out
}

.grsem_penalty <- function(x, lambda, alpha, spec) {
  group_part <- 0
  for (group_id in names(spec$group_weights)) {
    idx <- which(spec$groups == group_id & spec$penalty_factor > 0)
    if (length(idx)) {
      group_part <- group_part +
        spec$group_weights[[group_id]] * sqrt(sum(x[idx]^2))
    }
  }
  l1_idx <- which(spec$penalty_factor > 0 &
                    (!is.na(spec$groups) | spec$l1_only))
  l1_part <- sum(spec$l1_weights[l1_idx] * abs(x[l1_idx]))
  lambda * ((1 - alpha) * group_part + alpha * l1_part)
}

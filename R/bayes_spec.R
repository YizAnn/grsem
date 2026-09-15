#' Define Bayesian shrinkage settings
#'
#' @param family One of weak-information, regularized horseshoe, group
#'   horseshoe, or group-plus-local shrinkage.
#' @param global_scale Global shrinkage scale.
#' @param slab_scale Regularized-horseshoe slab scale.
#' @param rope Optional practical-equivalence interval.
#' @param sensitivity Named alternative prior settings.
#' @return A prior specification.
#' @export
grsem_prior <- function(family = c("weak", "regularized_horseshoe",
                                   "group_horseshoe", "group_local"),
                        global_scale = 0.1, slab_scale = 1,
                        rope = c(-0.05, 0.05), sensitivity = list()) {
  family <- match.arg(family)
  structure(list(
    family = family, global_scale = global_scale, slab_scale = slab_scale,
    rope = rope, sensitivity = sensitivity
  ), class = "grsem_prior_spec")
}

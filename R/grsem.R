#' Grouped regularized structural equation modeling
#'
#' Fits non-overlapping group-lasso and sparse-group-lasso paths while jointly
#' estimating all parameters in a continuous-variable, complete-data SEM.
#'
#' @param model lavaan model syntax.
#' @param data Complete numeric data frame.
#' @param groups Named parameter-to-group vector or selector data frame.
#' @param penalize Parameter classes eligible for regularization.
#' @param lambda Deprecated absolute lambda grid. Prefer `rho` so each sample
#'   computes its own `lambda_max`.
#' @param rho Relative path `lambda/lambda_max`; `NULL` creates an automatic
#'   path. Values must lie in `[0,1]`.
#' @param alpha Mixing grid: zero is group lasso, one is lasso.
#' @param penalty_factor Named or optimizer-length penalty factors. Zero means
#'   unpenalized. Grouped factors are zero/one in version 0.1.
#' @param group_weights Optional group weights; defaults to square root of group
#'   size.
#' @param estimator Only `"ML"` is supported in version 0.1.
#' @param std.lv Passed to lavaan.
#' @param warm_start Use the previous lambda solution as a starting value.
#' @param ungrouped Whether eligible but ungrouped parameters are unpenalized or
#'   treated as L1-only parameters.
#' @param meanstructure Passed to lavaan.
#' @param backend Estimation backend. `"lavaan_ml"` is the validated
#'   continuous backend; `"lavaan_dwls"` is an experimental ordered-indicator
#'   backend.
#' @param ordered Character vector of ordered indicators for the experimental
#'   categorical backend.
#' @param control Result of [grsem_control()].
#' @param ... Additional arguments passed to [lavaan::sem()].
#' @return Object of class `grsem_fit`.
#' @examples
#' # An overidentified CFA with four informative indicators and independent errors.
#' set.seed(1)
#' n <- 400
#' f <- rnorm(n)
#' d <- data.frame(
#'   x1 = 0.8 * f + rnorm(n, sd = 0.8),
#'   x2 = 0.9 * f + rnorm(n, sd = 0.8),
#'   x3 = 0.7 * f + rnorm(n, sd = 0.8),
#'   x4 = 0.85 * f + rnorm(n, sd = 0.8)
#' )
#' g <- c("f =~ x2"="candidate", "f =~ x3"="candidate",
#'        "f =~ x4"="candidate")
#' fit <- grsem("f =~ x1 + x2 + x3 + x4", d, groups=g,
#'              rho=c(0.2, 0), alpha=0,
#'              control=grsem_control(max_iter=100, tolerance=1e-3,
#'                                    n_starts=1))
#' fit
#' @export
grsem <- function(model, data, groups = NULL, penalize = c("loadings"),
                  lambda = NULL, rho = NULL, alpha = c(0, 0.25, 0.5),
                  penalty_factor = NULL, group_weights = NULL,
                  estimator = "ML", std.lv = TRUE, warm_start = TRUE,
                  ungrouped = c("unpenalized", "l1"),
                  meanstructure = FALSE,
                  backend = c("lavaan_ml", "lavaan_dwls"), ordered = NULL,
                  control = grsem_control(), ...) {
  ungrouped <- match.arg(ungrouped)
  backend <- match.arg(backend)
  if (!is.null(lambda) && !is.null(rho)) {
    stop("supply rho or the deprecated absolute lambda, not both", call. = FALSE)
  }
  if (!is.null(rho) && (any(!is.finite(rho)) || any(rho < 0 | rho > 1))) {
    stop("rho must be finite and lie in [0,1]", call. = FALSE)
  }
  allowed <- c("loadings", "regressions", "variances",
               "residual_covariances")
  if (!length(penalize) || any(!penalize %in% allowed)) {
    stop("penalize contains an unsupported parameter class", call. = FALSE)
  }
  if (any(!is.finite(alpha)) || any(alpha < 0 | alpha > 1)) {
    stop("alpha must lie in [0,1]", call. = FALSE)
  }
  backend_object <- .grsem_prepare_backend(
    model, data, estimator, std.lv, meanstructure,
    backend_type = backend, ordered = ordered, ...
  )
  spec <- .grsem_build_spec(
    backend_object, groups, penalize,
    penalty_factor, group_weights, ungrouped
  )
  if (!any(spec$penalty_factor > 0 &
           (!is.na(spec$groups) | spec$l1_only))) {
    stop("no parameters were selected for regularization", call. = FALSE)
  }

  restricted_x <- .grsem_restricted_start(backend_object, spec)
  automatic_path <- is.null(lambda) && is.null(rho)
  result_rows <- list()
  result_par <- list()
  row_id <- 0L

  for (a in unique(alpha)) {
    lambda_max <- .grsem_lambda_max(backend_object, spec, a, restricted_x)
    if (automatic_path) {
      rho_a <- if (lambda_max == 0) 0 else c(
        exp(seq(log(1), log(control$lambda_min_ratio),
                length.out = control$lambda_length)), 0
      )
      lambda_a <- rho_a * lambda_max
    } else if (!is.null(rho)) {
      rho_a <- sort(unique(as.numeric(rho)), decreasing = TRUE)
      lambda_a <- rho_a * lambda_max
    } else {
      sort(unique(as.numeric(lambda)), decreasing = TRUE)
      lambda_a <- sort(unique(as.numeric(lambda)), decreasing = TRUE)
      rho_a <- if (lambda_max > 0) lambda_a / lambda_max else
        ifelse(lambda_a == 0, 0, Inf)
    }
    if (any(!is.finite(lambda_a)) || any(lambda_a < 0)) {
      stop("lambda must be finite and non-negative", call. = FALSE)
    }
    current <- restricted_x
    for (path_index in seq_along(lambda_a)) {
      lam <- lambda_a[path_index]
      rho_value <- rho_a[path_index]
      row_id <- row_id + 1L
      if (lam == 0) {
        solution <- list(
          par = backend_object$x_unpenalized,
          raw_fit = backend_object$fn(backend_object$x_unpenalized),
          penalty = 0,
          objective = backend_object$fn(backend_object$x_unpenalized),
          status = "converged",
          converged = TRUE,
          approximately_converged = FALSE,
          model_valid = TRUE,
          prox_gradient_max = max(abs(backend_object$gr(
            backend_object$x_unpenalized))),
          relative_objective_change = 0,
          iterations = 0L,
          elapsed = 0,
          step = NA_real_,
          algorithm = "unpenalized",
          message = "unpenalized lavaan solution"
        )
      } else {
        starts <- list(if (warm_start) current else backend_object$x_unpenalized)
        if (control$n_starts >= 2L) starts[[2L]] <- backend_object$x_unpenalized
        if (control$n_starts >= 3L) {
          for (s in 3:control$n_starts) {
            starts[[s]] <- restricted_x +
              (s - 2) * 1e-3 * sign(backend_object$x_unpenalized - restricted_x)
          }
        }
        candidates <- lapply(seq_along(starts), function(start_id) {
          checkpoint_file <- NULL
          if (!is.null(control$checkpoint_path)) {
            checkpoint_file <- file.path(
              control$checkpoint_path,
              sprintf("alpha_%s_rho_%s_start_%02d.rds",
                      format(a, trim = TRUE),
                      format(rho_value, scientific = TRUE), start_id)
            )
          }
          .grsem_proximal_fit(
            starts[[start_id]], backend = backend_object, spec = spec,
            lambda = lam, alpha = a, control = control,
            checkpoint_file = checkpoint_file
          )
        })
        objectives <- vapply(candidates, `[[`, numeric(1), "objective")
        status_rank <- match(vapply(candidates, `[[`, character(1), "status"),
                             c("converged", "approximately_converged", "failed"))
        solution <- candidates[[order(status_rank, objectives)[1L]]]
      }
      current <- solution$par
      selected_parameter <- abs(solution$par) > control$zero_tolerance &
        spec$penalty_factor > 0
      selected_group <- vapply(names(spec$group_weights), function(group_id) {
        idx <- which(spec$groups == group_id & spec$penalty_factor > 0)
        any(abs(solution$par[idx]) > control$zero_tolerance)
      }, logical(1))
      k_eff <- sum(spec$penalty_factor == 0) + sum(selected_parameter)
      scale_fit <- 2 * backend_object$n * solution$raw_fit
      result_rows[[row_id]] <- data.frame(
        solution = row_id,
        alpha = a,
        rho = rho_value,
        lambda = lam,
        lambda_max = lambda_max,
        raw_fit = solution$raw_fit,
        penalty = solution$penalty,
        objective = solution$objective,
        status = solution$status,
        converged = solution$converged,
        approximately_converged = solution$approximately_converged,
        model_valid = solution$model_valid,
        prox_gradient_max = solution$prox_gradient_max,
        relative_objective_change = solution$relative_objective_change,
        iterations = solution$iterations,
        elapsed = solution$elapsed,
        algorithm = solution$algorithm,
        effective_parameters = k_eff,
        relative_aic_proxy = scale_fit + 2 * k_eff,
        relative_bic_proxy = scale_fit + log(backend_object$n) * k_eff,
        relative_ebic_proxy = scale_fit + log(backend_object$n) * k_eff +
          2 * 0.5 * log(max(1, length(solution$par))) * k_eff,
        nonzero_groups = sum(selected_group),
        nonzero_penalized = sum(selected_parameter),
        message = solution$message,
        stringsAsFactors = FALSE
      )
      result_par[[row_id]] <- solution$par
    }
  }

  object <- list(
    call = match.call(),
    model = model,
    data = data,
    backend = backend_object,
    backend_type = backend,
    spec = spec,
    path = do.call(rbind, result_rows),
    coefficients = result_par,
    control = control,
    limitations = c(
      if (backend == "lavaan_ml")
        "continuous variables, one group, complete data, ML only" else
        "experimental ordered-indicator DWLS backend; not a full-information categorical likelihood",
      "local stationary points only; global optimality is not claimed",
      "information criteria are explicitly labelled relative proxies",
      "ordinary chi-square, degrees of freedom, and standard errors are not assigned to penalized fits"
    ),
    refits = list()
  )
  class(object) <- "grsem_fit"
  object
}

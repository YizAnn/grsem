#' Optimizer controls for grsem
#'
#' @param algorithm Proximal-gradient algorithm. `"monotone_fista"` is the
#'   Phase 2 default; `"fista"` and `"pg"` are retained for comparison.
#' @param max_iter Maximum iterations.
#' @param tolerance Maximum proximal-gradient mapping tolerance required for
#'   strict convergence.
#' @param approximate_tolerance Maximum proximal-gradient mapping tolerance
#'   allowed for an approximately converged solution.
#' @param objective_tolerance Relative objective-change tolerance.
#' @param initial_step Initial gradient step size.
#' @param backtrack Multiplicative backtracking factor in `(0,1)`.
#' @param min_step Smallest permitted step.
#' @param max_backtrack Maximum backtracking attempts per iteration.
#' @param lambda_length Automatic relative-path length.
#' @param lambda_min_ratio Smallest automatic rho value.
#' @param zero_tolerance Threshold used for selection summaries.
#' @param n_starts Number of deterministic starting values per solution.
#' @param restart Restart FISTA when momentum is counterproductive.
#' @param checkpoint_path Optional directory for solution checkpoints.
#' @param checkpoint_every Save a checkpoint every this many iterations.
#' @param resume Resume compatible checkpoints when available.
#' @param keep_history Retain objective and first-order diagnostic histories.
#' @return Control list.
#' @export
grsem_control <- function(algorithm = c("monotone_fista", "fista", "pg"),
                          max_iter = 2000L, tolerance = 1e-5,
                          approximate_tolerance = 1e-4,
                          objective_tolerance = 1e-8,
                          initial_step = 1, backtrack = 0.5,
                          min_step = 1e-12, max_backtrack = 60L,
                          lambda_length = 30L, lambda_min_ratio = 0.01,
                          zero_tolerance = 1e-6, n_starts = 2L,
                          restart = TRUE, checkpoint_path = NULL,
                          checkpoint_every = 100L, resume = FALSE,
                          keep_history = FALSE) {
  algorithm <- match.arg(algorithm)
  if (!is.null(checkpoint_path) && (!is.character(checkpoint_path) ||
                                    length(checkpoint_path) != 1L)) {
    stop("checkpoint_path must be NULL or one directory", call. = FALSE)
  }
  list(
    algorithm = algorithm,
    max_iter = as.integer(max_iter),
    tolerance = tolerance,
    approximate_tolerance = approximate_tolerance,
    objective_tolerance = objective_tolerance,
    initial_step = initial_step,
    backtrack = backtrack,
    min_step = min_step,
    max_backtrack = as.integer(max_backtrack),
    lambda_length = as.integer(lambda_length),
    lambda_min_ratio = lambda_min_ratio,
    zero_tolerance = zero_tolerance,
    n_starts = as.integer(n_starts),
    restart = isTRUE(restart),
    checkpoint_path = checkpoint_path,
    checkpoint_every = as.integer(checkpoint_every),
    resume = isTRUE(resume),
    keep_history = isTRUE(keep_history)
  )
}

.grsem_restricted_start <- function(backend, spec) {
  penalized <- spec$penalty_factor > 0 &
    (!is.na(spec$groups) | spec$l1_only)
  x0 <- backend$x_unpenalized
  x0[penalized] <- 0
  free_u <- which(!penalized)
  if (!length(free_u)) return(x0)
  fn_u <- function(u) {
    x <- x0
    x[free_u] <- u
    backend$fn(x)
  }
  gr_u <- function(u) {
    x <- x0
    x[free_u] <- u
    backend$gr(x)[free_u]
  }
  opt <- stats::optim(
    par = x0[free_u], fn = fn_u, gr = gr_u,
    method = "L-BFGS-B",
    lower = backend$lower[free_u], upper = backend$upper[free_u],
    control = list(maxit = 2000, factr = 1e7)
  )
  x0[free_u] <- opt$par
  x0
}

.grsem_lambda_max <- function(backend, spec, alpha, restricted_x) {
  grad <- backend$gr(restricted_x)
  candidates <- 0
  group_ids <- names(spec$group_weights)
  for (group_id in group_ids) {
    idx <- which(spec$groups == group_id & spec$penalty_factor > 0)
    g <- grad[idx]
    w <- spec$group_weights[[group_id]]
    if (alpha == 0) {
      candidates <- c(candidates, sqrt(sum(g^2)) / w)
    } else if (alpha == 1) {
      candidates <- c(candidates,
                      max(abs(g) / spec$l1_weights[idx]))
    } else {
      condition <- function(lambda) {
        sqrt(sum(soft_threshold(
          g, lambda * alpha * spec$l1_weights[idx]
        )^2)) - lambda * (1 - alpha) * w
      }
      upper <- max(abs(g) / (alpha * spec$l1_weights[idx]))
      if (condition(0) > 0 && upper > 0) {
        candidates <- c(candidates, stats::uniroot(
          condition, interval = c(0, upper), tol = 1e-10
        )$root)
      }
    }
  }
  l1_idx <- which(spec$l1_only & spec$penalty_factor > 0)
  if (length(l1_idx)) {
    if (alpha == 0) {
      stop("ungrouped L1 parameters require alpha > 0", call. = FALSE)
    }
    candidates <- c(candidates,
                    abs(grad[l1_idx]) /
                      (alpha * spec$l1_weights[l1_idx]))
  }
  max(candidates[is.finite(candidates)], 0)
}

.grsem_prox_update <- function(x, gradient, step, lambda, alpha,
                               backend, spec) {
  candidate <- sparse_group_prox(
    x - step * gradient,
    groups = spec$groups,
    step_lambda = step * lambda,
    alpha = alpha,
    group_weights = spec$group_weights,
    penalty_factor = spec$l1_weights,
    l1_only = spec$l1_only
  )
  pmin(pmax(candidate, backend$lower), backend$upper)
}

.grsem_pg_norm <- function(x, step, backend, spec, lambda, alpha) {
  if (!is.finite(step) || step <= 0) return(Inf)
  grad <- backend$gr(x)
  if (any(!is.finite(grad))) return(Inf)
  prox <- .grsem_prox_update(x, grad, step, lambda, alpha, backend, spec)
  max(abs(x - prox)) / step
}

.grsem_line_search <- function(y, step, backend, spec, lambda, alpha, control) {
  raw_y <- backend$fn(y)
  grad_y <- backend$gr(y)
  if (!is.finite(raw_y) || any(!is.finite(grad_y))) {
    return(list(accepted = FALSE, message = "non-finite objective or gradient"))
  }
  for (bt in seq_len(control$max_backtrack)) {
    candidate <- .grsem_prox_update(
      y, grad_y, step, lambda, alpha, backend, spec
    )
    delta <- candidate - y
    raw_candidate <- backend$fn(candidate)
    majorizer <- raw_y + sum(grad_y * delta) + sum(delta^2) / (2 * step)
    if (is.finite(raw_candidate) && raw_candidate <= majorizer + 1e-12) {
      return(list(accepted = TRUE, par = candidate, raw = raw_candidate,
                  step = step, backtracks = bt - 1L))
    }
    step <- step * control$backtrack
    if (step < control$min_step) break
  }
  list(accepted = FALSE, message = "backtracking line search failed")
}

.grsem_checkpoint_signature <- function(lambda, alpha, dimension, algorithm) {
  paste(format(lambda, digits = 17), format(alpha, digits = 17),
        dimension, algorithm, sep = "|")
}

.grsem_save_checkpoint <- function(file, state) {
  if (is.null(file)) return(invisible(FALSE))
  file <- .grsem_output_path(file, create_parent = TRUE)
  saveRDS(state, file)
  invisible(TRUE)
}

.grsem_proximal_fit <- function(x0, backend, spec, lambda, alpha, control,
                                checkpoint_file = NULL) {
  if (!is.null(checkpoint_file)) checkpoint_file <- .grsem_output_path(checkpoint_file)
  start_time <- proc.time()[[3L]]
  algorithm <- control$algorithm %||% "pg"
  signature <- .grsem_checkpoint_signature(
    lambda, alpha, length(x0), algorithm
  )
  x <- pmin(pmax(x0, backend$lower), backend$upper)
  y <- x
  momentum <- 1
  step <- control$initial_step
  start_iteration <- 0L
  total <- backend$fn(x) + .grsem_penalty(x, lambda, alpha, spec)
  history <- list()
  if (isTRUE(control$resume) && !is.null(checkpoint_file) &&
      file.exists(checkpoint_file)) {
    saved <- try(readRDS(checkpoint_file), silent = TRUE)
    if (!inherits(saved, "try-error") && identical(saved$signature, signature)) {
      x <- saved$x
      y <- saved$y
      momentum <- saved$momentum
      step <- saved$step
      total <- saved$total
      start_iteration <- saved$iteration
      history <- saved$history
    }
  }

  message <- "maximum iterations reached"
  failure <- FALSE
  pg_max <- Inf
  rel_change <- Inf
  iteration <- start_iteration
  max_iteration <- control$max_iter

  if (start_iteration < max_iteration) {
    for (iteration in seq.int(start_iteration + 1L, max_iteration)) {
      search <- .grsem_line_search(
        y, step, backend, spec, lambda, alpha, control
      )
      if (!isTRUE(search$accepted)) {
        failure <- TRUE
        message <- search$message
        break
      }
      candidate <- search$par
      raw_candidate <- search$raw
      step <- search$step
      total_candidate <- raw_candidate +
        .grsem_penalty(candidate, lambda, alpha, spec)

      if (identical(algorithm, "monotone_fista") &&
          total_candidate > total + 1e-12) {
        y <- x
        momentum <- 1
        search <- .grsem_line_search(
          y, step, backend, spec, lambda, alpha, control
        )
        if (!isTRUE(search$accepted)) {
          failure <- TRUE
          message <- search$message
          break
        }
        candidate <- search$par
        raw_candidate <- search$raw
        step <- search$step
        total_candidate <- raw_candidate +
          .grsem_penalty(candidate, lambda, alpha, spec)
      }

      rel_change <- abs(total - total_candidate) / max(1, abs(total))
      pg_max <- .grsem_pg_norm(
        candidate, step, backend, spec, lambda, alpha
      )
      old_x <- x
      x <- candidate
      total <- total_candidate

      if (identical(algorithm, "pg")) {
        y <- x
        momentum <- 1
      } else {
        next_momentum <- (1 + sqrt(1 + 4 * momentum^2)) / 2
        y_candidate <- x + ((momentum - 1) / next_momentum) * (x - old_x)
        if (isTRUE(control$restart) && sum((x - old_x) *
                                          (y_candidate - x)) > 0) {
          y <- x
          momentum <- 1
        } else {
          y <- y_candidate
          momentum <- next_momentum
        }
      }

      if (isTRUE(control$keep_history)) {
        history[[length(history) + 1L]] <- c(
          iteration = iteration, objective = total,
          prox_gradient = pg_max, relative_change = rel_change,
          step = step
        )
      }
      if (!is.null(checkpoint_file) && control$checkpoint_every > 0L &&
          iteration %% control$checkpoint_every == 0L) {
        .grsem_save_checkpoint(checkpoint_file, list(
          signature = signature, x = x, y = y, momentum = momentum,
          step = step, total = total, iteration = iteration,
          history = history
        ))
      }
      if (pg_max <= control$tolerance &&
          rel_change <= control$objective_tolerance &&
          iteration < max_iteration) {
        message <- "strict stopping rules reached"
        break
      }
      step <- min(step / control$backtrack, control$initial_step)
    }
  }

  raw <- backend$fn(x)
  penalty <- .grsem_penalty(x, lambda, alpha, spec)
  total <- raw + penalty
  model_valid <- is.finite(raw) && is.finite(total) &&
    all(is.finite(x)) && (is.null(backend$valid) || isTRUE(backend$valid(x)))
  hit_limit <- iteration >= max_iteration
  strict <- model_valid && !failure && !hit_limit &&
    pg_max <= control$tolerance &&
    rel_change <= control$objective_tolerance
  approximate <- model_valid && !failure && !strict &&
    pg_max <= control$approximate_tolerance
  status <- if (strict) {
    "converged"
  } else if (approximate) {
    "approximately_converged"
  } else {
    "failed"
  }
  if (status == "approximately_converged" && hit_limit) {
    message <- "iteration limit reached with small first-order error"
  } else if (status == "failed" && hit_limit) {
    message <- "iteration limit reached without acceptable first-order error"
  }
  .grsem_save_checkpoint(checkpoint_file, list(
    signature = signature, x = x, y = y, momentum = momentum,
    step = step, total = total, iteration = iteration,
    history = history, complete = status != "failed", status = status
  ))

  list(
    par = x,
    raw_fit = raw,
    penalty = penalty,
    objective = total,
    status = status,
    converged = identical(status, "converged"),
    approximately_converged = identical(status, "approximately_converged"),
    model_valid = model_valid,
    prox_gradient_max = pg_max,
    relative_objective_change = rel_change,
    iterations = iteration,
    elapsed = proc.time()[[3L]] - start_time,
    step = step,
    algorithm = algorithm,
    message = message,
    history = if (isTRUE(control$keep_history)) do.call(rbind, history) else NULL
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

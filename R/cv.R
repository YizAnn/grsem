.grsem_gaussian_nll <- function(sigma, mu, newdata) {
  sigma <- as.matrix(sigma)
  vars <- colnames(sigma)
  if (!is.numeric(sigma) || !nrow(sigma) || nrow(sigma) != ncol(sigma) ||
      is.null(vars) || anyNA(vars) || anyDuplicated(vars) ||
      !identical(rownames(sigma), vars) || any(!is.finite(sigma)) ||
      !isTRUE(isSymmetric(sigma))) stop("invalid named covariance matrix")
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (!nrow(newdata) || anyDuplicated(names(newdata)) ||
      !all(vars %in% names(newdata))) stop("test data must contain all model variables")
  if (!all(vapply(newdata[vars], is.numeric, logical(1)))) {
    stop("test variables must be numeric")
  }
  y <- as.matrix(newdata[, vars, drop = FALSE])
  if (any(!is.finite(y))) stop("test data contain missing or non-finite values")
  if (is.null(mu)) {
    mu <- stats::setNames(rep(0, length(vars)), vars)
  } else if (is.null(names(mu))) {
    if (length(mu) != length(vars)) stop("unnamed mean has wrong length")
    names(mu) <- vars
  } else {
    if (anyDuplicated(names(mu)) || !all(vars %in% names(mu))) {
      stop("mean must name all model variables exactly once")
    }
    mu <- mu[vars]
  }
  if (!is.numeric(mu) || any(!is.finite(mu))) stop("invalid model mean")
  centered <- sweep(y, 2L, mu[vars], "-")
  chol_sigma <- chol(sigma)
  whitened <- backsolve(chol_sigma, t(centered), transpose = TRUE)
  mean(0.5 * (length(vars) * log(2 * pi) +
                2 * sum(log(diag(chol_sigma))) + colSums(whitened^2)))
}

.grsem_test_nll <- function(fit, newdata) {
  implied <- lavaan::lavInspect(fit, "implied")
  .grsem_gaussian_nll(implied$cov, implied$mean, newdata)
}

#' K-fold evaluation with explicit failures
#' @inheritParams grsem
#' @param folds Number of folds, between two and the number of rows.
#' @param seed Random seed.
#' @return One row per fold and candidate, including evaluation_status,
#'   failure_stage and message. No failed folds are dropped or averaged away.
#' @export
cv_grsem <- function(model, data, groups = NULL,
                     penalize = "loadings", rho = NULL,
                     alpha = c(0, 0.5, 1), penalty_factor = NULL,
                     group_weights = NULL, estimator = "ML", std.lv = TRUE,
                     ungrouped = "unpenalized", meanstructure = FALSE,
                     control = grsem_control(), folds = 5L, seed = 1L, ...) {
  if (length(folds) != 1L || !is.finite(folds) || folds != as.integer(folds) ||
      folds < 2 || folds > nrow(data)) stop("invalid folds", call. = FALSE)
  if (!length(alpha) || any(!is.finite(alpha)) || any(alpha < 0 | alpha > 1)) {
    stop("invalid alpha", call. = FALSE)
  }
  if (is.null(rho)) rho <- c(exp(seq(log(1), log(control$lambda_min_ratio),
                                    length.out = control$lambda_length)), 0)
  if (!length(rho) || any(!is.finite(rho)) || any(rho < 0 | rho > 1)) {
    stop("invalid rho", call. = FALSE)
  }
  rho <- sort(unique(rho), decreasing = TRUE)
  alpha <- unique(alpha)
  grid <- expand.grid(rho = rho, alpha = alpha)
  set.seed(seed)
  fold_id <- sample(rep(seq_len(folds), length.out = nrow(data)))
  out <- list()
  for (fold in seq_len(folds)) {
    train <- data[fold_id != fold, , drop = FALSE]
    test <- data[fold_id == fold, , drop = FALSE]
    fit <- tryCatch(grsem(
      model = model, data = train, groups = groups, penalize = penalize,
      rho = rho, alpha = alpha, penalty_factor = penalty_factor,
      group_weights = group_weights, estimator = estimator, std.lv = std.lv,
      warm_start = TRUE, ungrouped = ungrouped, meanstructure = meanstructure,
      control = control, ...
    ), error = function(e) e)
    for (j in seq_len(nrow(grid))) {
      row <- data.frame(
        fold = fold, solution = j, alpha = grid$alpha[j], rho = grid$rho[j],
        lambda_max = NA_real_, lambda = NA_real_, status = "failed",
        nonzero_groups = NA_integer_, nonzero_penalized = NA_integer_,
        prox_gradient_max = NA_real_, objective = NA_real_, test_nll = NA_real_,
        evaluation_status = "failed", failure_stage = "fit", message = "",
        stringsAsFactors = FALSE
      )
      if (inherits(fit, "error")) {
        row$message <- conditionMessage(fit)
      } else {
        k <- which(fit$path$alpha == row$alpha & fit$path$rho == row$rho)
        if (length(k) != 1L) {
          row$message <- "candidate missing or duplicated in fitted path"
        } else {
          for (nm in c("solution", "lambda_max", "lambda", "status",
                       "nonzero_groups", "nonzero_penalized",
                       "prox_gradient_max", "objective")) row[[nm]] <- fit$path[[nm]][k]
          if (row$status == "failed") {
            row$message <- fit$path$message[k]
          } else {
            stage <- "refit"
            result <- tryCatch({
              refit <- post_selection_refit(fit, row$solution, data = train)
              if (!isTRUE(lavaan::lavInspect(refit, "converged"))) {
                stop("post-selection refit did not converge")
              }
              stage <- "nll"
              nll <- .grsem_test_nll(refit, test)
              if (!is.finite(nll)) stop("non-finite test NLL")
              nll
            }, error = function(e) e)
            if (inherits(result, "error")) {
              row$failure_stage <- stage
              row$message <- conditionMessage(result)
            } else {
              row$test_nll <- result
              row$evaluation_status <- "ok"
              row$failure_stage <- ""
            }
          }
        }
      }
      out[[length(out) + 1L]] <- row
    }
  }
  do.call(rbind, out)
}

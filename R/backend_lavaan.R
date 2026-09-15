.grsem_check_backend <- function() {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("lavaan is required", call. = FALSE)
  }
  if (utils::packageVersion("lavaan") < "0.7.2") {
    stop("grsem requires lavaan >= 0.7-2 for lav_export_estimation()",
         call. = FALSE)
  }
}

.grsem_prepare_backend <- function(model, data, estimator, std.lv,
                                   meanstructure,
                                   backend_type = c("lavaan_ml", "lavaan_dwls"),
                                   ordered = NULL, ...) {
  .grsem_check_backend()
  backend_type <- match.arg(backend_type)
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (anyNA(data)) {
    stop("version 0.1 supports complete data only; FIML is a later milestone",
         call. = FALSE)
  }
  if (backend_type == "lavaan_ml") {
    if (any(!vapply(data, is.numeric, logical(1)))) {
      stop("lavaan_ml supports numeric variables only", call. = FALSE)
    }
    if (!identical(toupper(estimator), "ML")) {
      stop("lavaan_ml supports estimator = 'ML' only", call. = FALSE)
    }
    fit <- lavaan::sem(
      model = model, data = data, estimator = "ML", std.lv = std.lv,
      meanstructure = meanstructure, se = "none", test = "none", ...
    )
  } else {
    if (is.null(ordered) || !length(ordered) ||
        any(!ordered %in% names(data))) {
      stop("lavaan_dwls requires a valid ordered variable vector", call. = FALSE)
    }
    if (any(vapply(data[ordered], function(x) length(unique(x)) > 2L,
                   logical(1)))) {
      stop("the experimental backend currently supports binary indicators only",
           call. = FALSE)
    }
    fit <- lavaan::sem(
      model = model, data = data, estimator = "DWLS", ordered = ordered,
      parameterization = "theta", std.lv = std.lv,
      meanstructure = TRUE, se = "none", test = "none", ...
    )
  }
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop("the unpenalized lavaan template did not converge", call. = FALSE)
  }
  exported <- lavaan::lav_export_estimation(fit)
  n_native <- length(exported$starting_values)
  pt <- lavaan::parameterTable(fit)
  free_ids <- sort(unique(pt$free[pt$free > 0L]))
  n_full <- length(free_ids)
  to_native_matrix <- diag(n_native)
  full_to_opt <- seq_len(n_native)
  full_scale <- rep(1, n_native)
  if (fit@Model@eq.constraints) {
    K <- fit@Model@eq.constraints.K
    if (nrow(K) != n_full || ncol(K) != n_native) {
      stop("unsupported lavaan equality-constraint dimensions", call. = FALSE)
    }
    normalized <- matrix(0, nrow(K), ncol(K))
    scales <- numeric(nrow(K))
    for (i in seq_len(nrow(K))) {
      ni <- sqrt(sum(K[i, ]^2))
      if (!is.finite(ni) || ni <= 1e-12) {
        stop("unsupported zero row in equality-constraint map", call. = FALSE)
      }
      u <- K[i, ] / ni
      first <- which(abs(u) > 1e-10)[1]
      if (u[first] < 0) u <- -u
      normalized[i, ] <- u
      scales[i] <- sum(K[i, ] * u)
    }
    key <- apply(round(normalized, 10), 1, paste, collapse = ":")
    keep <- !duplicated(key)
    U <- normalized[keep, , drop = FALSE]
    if (nrow(U) != n_native ||
        max(abs(U %*% t(U) - diag(n_native))) > 1e-7) {
      stop("version 0.1 supports only simple orthogonal equality constraints",
           call. = FALSE)
    }
    full_to_opt <- match(key, key[keep])
    full_scale <- scales
    to_native_matrix <- t(U)
  }
  to_native <- function(x) as.numeric(to_native_matrix %*% x)
  from_native <- function(x) as.numeric(t(to_native_matrix) %*% x)
  npar <- n_native
  lower <- exported$lower
  upper <- exported$upper
  if (is.null(lower)) lower <- rep(-Inf, npar)
  if (is.null(upper)) upper <- rep( Inf, npar)
  lower <- rep_len(lower, npar)
  upper <- rep_len(upper, npar)
  if (fit@Model@eq.constraints &&
      (any(is.finite(lower)) || any(is.finite(upper)))) {
    stop("rotated equality constraints with finite parameter bounds are not yet supported",
         call. = FALSE)
  }
  objective <- function(x) exported$objective_function(to_native(x), fit)
  gradient <- function(x) as.numeric(t(to_native_matrix) %*%
      exported$gradient_function(to_native(x), fit))
  list(
    fit = fit,
    exported = exported,
    fn = objective,
    gr = gradient,
    valid = function(x) {
      value <- try(objective(x), silent = TRUE)
      !inherits(value, "try-error") && length(value) == 1L && is.finite(value)
    },
    get_coef = function(x) exported$get_coef(to_native(x), fit),
    lower = lower,
    upper = upper,
    x_unpenalized = from_native(exported$starting_values),
    n = nrow(data),
    data = data,
    model = model,
    std.lv = std.lv,
    meanstructure = meanstructure,
    estimator = if (backend_type == "lavaan_ml") "ML" else "DWLS",
    backend_type = backend_type,
    ordered = ordered,
    full_to_opt = full_to_opt,
    full_scale = full_scale,
    to_native = to_native
  )
}

.grsem_parameter_table <- function(object, x) {
  pt <- object$spec$full_parameter_table
  coef_full <- if (is.function(object$backend$get_coef)) {
    object$backend$get_coef(x)
  } else {
    object$backend$exported$get_coef(x, object$backend$fit)
  }
  free_ids <- sort(unique(pt$free[pt$free > 0L]))
  est_by_free <- rep(NA_real_, max(free_ids))
  est_by_free[free_ids] <- as.numeric(coef_full)
  pt$est <- ifelse(pt$free > 0L, est_by_free[pmax(pt$free, 1L)], pt$ustart)
  pt[, intersect(c("lhs", "op", "rhs", "label", "free", "est"), names(pt)),
     drop = FALSE]
}

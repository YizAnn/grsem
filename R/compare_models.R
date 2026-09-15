.grsem_comparison_row <- function(name) {
  data.frame(
    model = name, status = "failed", CFI = NA_real_, TLI = NA_real_,
    RMSEA = NA_real_, SRMR = NA_real_, negative_variances = NA_integer_,
    standardized_paths_over_one = NA_integer_,
    test_common_covariance_rmse = NA_real_,
    ordinary_information_criteria_comparable = FALSE, message = "",
    stringsAsFactors = FALSE
  )
}

#' Compare fitted models with a fixed result schema
#' @param fits Named list of lavaan fits or post-selection refits.
#' @param test_data Optional held-out data frame or named list.
#' @param common_variables Exact common observed variables for all comparisons.
#' @return A data frame, preserving failed models and error messages.
#' @export
compare_models <- function(fits, test_data = NULL, common_variables = NULL) {
  if (!is.list(fits)) stop("fits must be a list")
  if (!length(fits)) return(.grsem_comparison_row("")[FALSE, ])
  if (is.null(names(fits))) names(fits) <- paste0("model", seq_along(fits))
  if (anyNA(names(fits)) || any(!nzchar(names(fits))) || anyDuplicated(names(fits))) {
    stop("fits must have unique nonempty names")
  }
  rows <- lapply(seq_along(fits), function(i) {
    row <- .grsem_comparison_row(names(fits)[i])
    tryCatch({
      fit <- fits[[i]]
      if (!inherits(fit, "lavaan")) stop("not a fitted lavaan post-selection model")
      if (!isTRUE(lavaan::lavInspect(fit, "converged"))) stop("model did not converge")
      fm <- lavaan::fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
      row[1, c("CFI", "TLI", "RMSEA", "SRMR")] <- as.list(unname(fm))
      if (!is.null(test_data) && !is.null(common_variables)) {
        data_now <- if (is.list(test_data) && !is.data.frame(test_data)) {
          test_data[[names(fits)[i]]]
        } else test_data
        sigma <- lavaan::lavInspect(fit, "implied")$cov
        vars <- common_variables
        if (length(vars) < 2L || anyNA(vars) || anyDuplicated(vars) ||
            !all(vars %in% colnames(sigma)) || !all(vars %in% names(data_now))) {
          stop("all common_variables must be present in model and test data")
        }
        if (nrow(data_now) < 2L ||
            !all(vapply(data_now[vars], is.numeric, logical(1))) ||
            any(!is.finite(as.matrix(data_now[vars])))) stop("invalid test data")
        row$test_common_covariance_rmse <- sqrt(mean(
          (stats::cov(data_now[, vars, drop = FALSE]) - sigma[vars, vars, drop = FALSE])^2
        ))
      }
      pe <- lavaan::standardizedSolution(fit)
      pv <- lavaan::parameterEstimates(fit)
      row$negative_variances <- sum(pv$op == "~~" & pv$lhs == pv$rhs & pv$est < 0)
      row$standardized_paths_over_one <- sum(pe$op == "~" & abs(pe$est.std) > 1)
      row$status <- "converged"
      row$message <- "Information criteria are not mechanically compared across item sets."
      row
    }, error = function(e) {
      failed <- .grsem_comparison_row(names(fits)[i])
      failed$message <- conditionMessage(e)
      failed
    })
  })
  do.call(rbind, rows)
}

#' Define a unified grsem model specification
#'
#' The specification is shared by penalized estimation, candidate-role
#' diagnostics, and experimental Bayesian backends. It stores definitions, not
#' fitted factor scores.
#'
#' @param model lavaan model syntax.
#' @param data Analysis data frame.
#' @param module_map Named item-to-module vector or data frame with `item` and
#'   `module` columns.
#' @param ordered Optional ordered/binary indicator names.
#' @param constraints Optional theory constraints for role learning.
#' @param metadata Optional provenance and version metadata.
#' @return A `grsem_model_spec` object.
#' @export
grsem_model <- function(model, data, module_map = NULL, ordered = NULL,
                        constraints = list(), metadata = list()) {
  if (!is.character(model) || length(model) != 1L) {
    stop("model must be one lavaan syntax string", call. = FALSE)
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (!is.null(ordered) && any(!ordered %in% names(data))) {
    stop("ordered variables were not found in data", call. = FALSE)
  }
  if (is.data.frame(module_map)) {
    if (!all(c("item", "module") %in% names(module_map))) {
      stop("module_map data frame needs item and module columns", call. = FALSE)
    }
    module_map <- stats::setNames(as.character(module_map$module),
                                  as.character(module_map$item))
  }
  if (!is.null(module_map)) {
    if (is.null(names(module_map)) || any(!names(module_map) %in% names(data))) {
      stop("module_map must be a named item-to-module vector", call. = FALSE)
    }
    module_map <- stats::setNames(as.character(module_map), names(module_map))
  }
  object <- list(
    model = model,
    data = data,
    module_map = module_map,
    ordered = ordered,
    constraints = constraints,
    metadata = metadata,
    created = Sys.time(),
    schema_version = "0.2"
  )
  class(object) <- "grsem_model_spec"
  object
}

#' @export
print.grsem_model_spec <- function(x, ...) {
  cat("grsem model specification\n")
  cat("Rows:", nrow(x$data), " Variables:", ncol(x$data), "\n")
  cat("Ordered indicators:", length(x$ordered), "\n")
  cat("Business modules:", length(unique(x$module_map)), "\n")
  invisible(x)
}

#' Export selected structural regression edges
#'
#' Only op == "~" is a directed edge (rhs -> lhs). Nonzero fixed and
#' unpenalized regressions are included and identified in metadata. Loadings,
#' variances and covariances are excluded. This is not causal discovery.
#' @param object A grsem_fit.
#' @param solution Explicit solution identifier; failed solutions are rejected.
#' @param zero_tolerance Nonnegative threshold; defaults to the fit control.
#'   Only abs(estimate) > zero_tolerance is included, in exported parameter units.
#' @return Data frame with from, to, estimate and parameter/solution metadata.
#' @export
selected_edges <- function(object, solution, zero_tolerance = NULL) {
  if (!inherits(object, "grsem_fit")) stop("object must be a grsem_fit")
  if (length(solution) != 1L || is.na(solution)) stop("supply one solution")
  row <- match(solution, object$path$solution)
  if (is.na(row)) stop("unknown solution")
  if (object$path$status[row] == "failed") stop("cannot export a failed solution")
  if (is.null(zero_tolerance)) zero_tolerance <- object$control$zero_tolerance
  if (length(zero_tolerance) != 1L || !is.finite(zero_tolerance) ||
      zero_tolerance < 0) stop("invalid zero_tolerance")
  pt <- .grsem_parameter_table(object, object$coefficients[[row]])
  pt$parameter_id <- seq_len(nrow(pt))
  free_ids <- sort(unique(object$spec$full_parameter_table$free[
    object$spec$full_parameter_table$free > 0L]))
  opt <- object$spec$full_to_opt[match(pt$free, free_ids)]
  penalized <- rep(FALSE, nrow(pt))
  group <- rep(NA_character_, nrow(pt))
  mapped <- which(!is.na(opt) & opt > 0L)
  penalized[mapped] <- object$spec$penalty_factor[opt[mapped]] > 0
  group[mapped] <- object$spec$groups[opt[mapped]]
  structural <- pt$op == "~"
  if (any(!is.finite(pt$est[structural]))) stop("non-finite structural estimate")
  keep <- which(structural & abs(pt$est) > zero_tolerance)
  data.frame(
    from = pt$rhs[keep], to = pt$lhs[keep], estimate = pt$est[keep],
    parameter_id = pt$parameter_id[keep], lhs = pt$lhs[keep], op = pt$op[keep],
    rhs = pt$rhs[keep], label = pt$label[keep], free = pt$free[keep],
    is_fixed = pt$free[keep] == 0L, penalized = penalized[keep], group = group[keep],
    solution = rep(solution, length(keep)),
    status = rep(object$path$status[row], length(keep)),
    zero_tolerance = rep(zero_tolerance, length(keep)),
    stringsAsFactors = FALSE
  )
}

#' Diagnose acyclicity of a selected structural regression graph
#'
#' This checks graph legality only, not identification or causal validity.
#' @param edges Data frame with nonempty character from/to columns.
#' @param nodes Optional additional nodes, including isolated vertices.
#' @return List with is_dag, topological_order (DAG only), cycle (closed
#'   directed walk for a non-DAG), and nodes. Empty graphs are DAGs.
#' @export
validate_dag <- function(edges, nodes = character()) {
  if (!is.data.frame(edges) || !all(c("from", "to") %in% names(edges))) {
    stop("edges must be a data frame with from and to")
  }
  valid_names <- function(x) is.character(x) && !anyNA(x) && all(nzchar(x))
  if (!valid_names(edges$from) || !valid_names(edges$to) || !valid_names(nodes)) {
    stop("node names must be nonempty, nonmissing character values")
  }
  nodes <- unique(c(nodes, edges$from, edges$to))
  pairs <- unique(edges[, c("from", "to"), drop = FALSE])
  adjacency <- lapply(nodes, function(x) pairs$to[pairs$from == x])
  names(adjacency) <- nodes
  indegree <- stats::setNames(integer(length(nodes)), nodes)
  for (x in pairs$to) indegree[x] <- indegree[x] + 1L
  queue <- nodes[indegree == 0L]
  order <- character()
  while (length(queue)) {
    x <- queue[1L]
    queue <- queue[-1L]
    order <- c(order, x)
    for (y in adjacency[[x]]) {
      indegree[y] <- indegree[y] - 1L
      if (indegree[y] == 0L) queue <- c(queue, y)
    }
  }
  is_dag <- length(order) == length(nodes)
  cycle <- character()
  if (!is_dag) {
    # Iterative DFS avoids a recursion-depth limit for long graphs.
    color <- stats::setNames(integer(length(nodes)), nodes)
    for (start in nodes) {
      if (color[start] != 0L || length(cycle)) next
      stack <- start
      next_child <- 1L
      color[start] <- 1L
      while (length(stack) && !length(cycle)) {
        depth <- length(stack)
        x <- stack[depth]
        children <- adjacency[[x]]
        k <- next_child[depth]
        if (k > length(children)) {
          color[x] <- 2L
          stack <- stack[-depth]
          next_child <- next_child[-depth]
        } else {
          next_child[depth] <- k + 1L
          y <- children[k]
          if (color[y] == 1L) {
            cycle <- c(stack[seq.int(match(y, stack), depth)], y)
          } else if (color[y] == 0L) {
            color[y] <- 1L
            stack <- c(stack, y)
            next_child <- c(next_child, 1L)
          }
        }
      }
    }
  }
  list(is_dag = is_dag, topological_order = if (is_dag) order else character(),
       cycle = cycle, nodes = nodes)
}

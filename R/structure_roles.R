.grsem_role_data <- function(spec, selected_items, module_map, constraints) {
  if (is.null(module_map)) stop("module_map is required", call. = FALSE)
  selected_items <- intersect(selected_items, names(module_map))
  if (!length(selected_items)) stop("no selected items match module_map", call. = FALSE)
  modules <- split(selected_items, module_map[selected_items])
  out <- lapply(modules, function(items) {
    rowMeans(as.data.frame(lapply(spec$data[items], as.numeric)), na.rm = FALSE)
  })
  names(out) <- paste0("Module_", make.names(names(out), unique = TRUE))
  core <- constraints$core_nodes
  if (!is.null(core)) {
    if (!is.list(core) || is.null(names(core))) {
      stop("constraints$core_nodes must be a named item-list", call. = FALSE)
    }
    for (node in names(core)) {
      items <- core[[node]]
      if (any(!items %in% names(spec$data))) {
        stop("core-node items were not found: ", node, call. = FALSE)
      }
      out[[node]] <- rowMeans(
        as.data.frame(lapply(spec$data[items], as.numeric)), na.rm = FALSE
      )
    }
  }
  data <- as.data.frame(out, check.names = FALSE)
  data <- data[stats::complete.cases(data), , drop = FALSE]
  data[] <- lapply(data, function(x) as.numeric(scale(x)))
  data
}

.grsem_fixed_gaps <- function(labels, constraints) {
  gaps <- matrix(FALSE, length(labels), length(labels),
                 dimnames = list(labels, labels))
  pairs <- constraints$forbidden_adjacent
  if (!is.null(pairs)) {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (ncol(pairs) < 2L) stop("forbidden_adjacent needs two columns", call. = FALSE)
    for (i in seq_len(nrow(pairs))) {
      a <- as.character(pairs[i, 1]); b <- as.character(pairs[i, 2])
      if (a %in% labels && b %in% labels) gaps[a, b] <- gaps[b, a] <- TRUE
    }
  }
  gaps
}

.grsem_matrix_edges <- function(amat, labels, algorithm, replicate) {
  amat <- as.matrix(amat)
  dimnames(amat) <- list(labels, labels)
  rows <- list()
  if (length(labels) < 2L) return(data.frame())
  for (i in seq_len(length(labels) - 1L)) for (j in (i + 1L):length(labels)) {
    a <- amat[i, j]; b <- amat[j, i]
    if (a == 0 && b == 0) next
    direction <- if (algorithm == "fci") {
      if (a == 2 && b == 1) paste(labels[i], "->", labels[j]) else
        if (a == 1 && b == 2) paste(labels[j], "->", labels[i]) else
          "unidentified"
    } else if (a == 1 && b == 0) {
      paste(labels[i], "->", labels[j])
    } else if (a == 0 && b == 1) {
      paste(labels[j], "->", labels[i])
    } else {
      "unidentified"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      algorithm = algorithm, replicate = replicate,
      node1 = labels[i], node2 = labels[j],
      skeleton = paste(sort(c(labels[i], labels[j])), collapse = " -- "),
      direction = direction, direction_identified = direction != "unidentified",
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

.grsem_run_structure_method <- function(data, method, constraints,
                                         alpha, m_max, replicate) {
  labels <- names(data)
  if (method == "notears") {
    return(list(status = "not_implemented", edges = data.frame(),
                message = "NOTEARS backend is not installed in the Phase 2 environment"))
  }
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    return(list(status = "not_implemented", edges = data.frame(),
                message = "pcalg is required"))
  }
  gaps <- .grsem_fixed_gaps(labels, constraints)
  suff <- list(C = stats::cor(data), n = nrow(data))
  result <- try({
    if (method == "pc") {
      fit <- pcalg::pc(
        suff, pcalg::gaussCItest, alpha = alpha, labels = labels,
        fixedGaps = gaps, skel.method = "stable", m.max = m_max,
        u2pd = "relaxed", solve.confl = TRUE, verbose = FALSE
      )
      .grsem_matrix_edges(as(fit@graph, "matrix"), labels, method, replicate)
    } else if (method == "fci") {
      fit <- pcalg::fci(
        suff, pcalg::gaussCItest, alpha = alpha, labels = labels,
        fixedGaps = gaps, skel.method = "stable", m.max = m_max,
        pdsep.max = m_max, verbose = FALSE
      )
      .grsem_matrix_edges(fit@amat, labels, method, replicate)
    } else if (method == "ges") {
      score <- methods::new(
        "GaussL0penObsScore", data = as.matrix(data),
        lambda = 0.5 * log(nrow(data)), intercept = FALSE
      )
      fit <- pcalg::ges(
        score, labels = labels, fixedGaps = gaps,
        phase = c("forward", "backward"), iterate = FALSE, verbose = FALSE
      )
      .grsem_matrix_edges(as(fit$essgraph, "matrix"), labels, method, replicate)
    } else {
      stop("unknown method: ", method)
    }
  }, silent = TRUE)
  if (inherits(result, "try-error")) {
    list(status = "failed", edges = data.frame(), message = as.character(result))
  } else {
    list(status = "converged", edges = result, message = "")
  }
}

.grsem_theory_direction <- function(node1, node2, constraints) {
  tiers <- constraints$tiers
  if (is.null(tiers) || is.null(names(tiers)) ||
      !all(c(node1, node2) %in% names(tiers)) || tiers[node1] == tiers[node2]) {
    return("unidentified")
  }
  if (tiers[node1] < tiers[node2]) paste(node1, "->", node2) else
    paste(node2, "->", node1)
}

#' Learn candidate structural roles under theory restrictions
#'
#' Regularization supplies candidate items first. This auxiliary routine then
#' learns module-level skeletons; it does not assign business value or claim a
#' unique causal graph.
#'
#' @param model A [grsem_model()] specification.
#' @param selected_items Training-selected topic items.
#' @param module_map Optional item-to-module map; defaults to the model spec.
#' @param constraints Theory tiers, core-node definitions, and optional
#'   forbidden adjacencies.
#' @param methods Any of `pc`, `fci`, `ges`, and `notears`.
#' @param bootstrap Number of respondent bootstrap resamples.
#' @param alpha Conditional-independence test level.
#' @param m_max Maximum conditioning-set size.
#' @param seed Random seed.
#' @return A `grsem_role_learning` object.
#' @export
learn_candidate_roles <- function(model, selected_items,
                                  module_map = model$module_map,
                                  constraints = model$constraints,
                                  methods = c("pc", "fci", "ges", "notears"),
                                  bootstrap = 100L, alpha = 0.05,
                                  m_max = 3L, seed = 1L) {
  if (!inherits(model, "grsem_model_spec")) {
    stop("model must be created by grsem_model()", call. = FALSE)
  }
  data <- .grsem_role_data(model, selected_items, module_map, constraints)
  methods <- unique(tolower(methods))
  set.seed(seed)
  edge_rows <- list(); task_rows <- list()
  for (method in methods) {
    for (b in 0:as.integer(bootstrap)) {
      sample_data <- if (b == 0L) data else
        data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
      result <- .grsem_run_structure_method(
        sample_data, method, constraints, alpha, m_max, b
      )
      task_rows[[length(task_rows) + 1L]] <- data.frame(
        algorithm = method, replicate = b, status = result$status,
        edge_count = nrow(result$edges), message = result$message,
        stringsAsFactors = FALSE
      )
      if (nrow(result$edges)) edge_rows[[length(edge_rows) + 1L]] <- result$edges
    }
  }
  edges <- if (length(edge_rows)) do.call(rbind, edge_rows) else data.frame()
  tasks <- do.call(rbind, task_rows)
  boot <- edges[edges$replicate > 0L, , drop = FALSE]
  frequencies <- if (!nrow(boot)) data.frame() else {
    keys <- interaction(boot$algorithm, boot$skeleton, drop = TRUE)
    do.call(rbind, lapply(split(boot, keys), function(z) {
      successful <- sum(tasks$algorithm == z$algorithm[1] &
                          tasks$replicate > 0L & tasks$status == "converged")
      data.frame(
        algorithm = z$algorithm[1], node1 = z$node1[1], node2 = z$node2[1],
        skeleton = z$skeleton[1], successful_bootstraps = successful,
        adjacency_frequency = length(unique(z$replicate)) / max(1, successful),
        identified_direction_frequency = mean(z$direction_identified),
        theory_compatible_direction = .grsem_theory_direction(
          z$node1[1], z$node2[1], constraints
        ), stringsAsFactors = FALSE
      )
    }))
  }
  main <- edges[edges$replicate == 0L, , drop = FALSE]
  algorithms <- unique(main$algorithm)
  agreement <- list()
  if (length(algorithms) >= 2L) {
    for (i in seq_len(length(algorithms) - 1L)) for (j in (i + 1L):length(algorithms)) {
      a <- unique(main$skeleton[main$algorithm == algorithms[i]])
      b <- unique(main$skeleton[main$algorithm == algorithms[j]])
      agreement[[length(agreement) + 1L]] <- data.frame(
        algorithm1 = algorithms[i], algorithm2 = algorithms[j],
        skeleton_jaccard = length(intersect(a, b)) / max(1, length(union(a, b)))
      )
    }
  }
  structure(list(
    data = data, main_edges = main, bootstrap_edges = boot,
    frequencies = frequencies, tasks = tasks,
    agreement = if (length(agreement)) do.call(rbind, agreement) else data.frame(),
    constraints = constraints,
    interpretation = "Skeleton and stable adjacency are primary; unidentified directions remain unidentified."
  ), class = "grsem_role_learning")
}

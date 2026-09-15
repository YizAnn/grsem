.grsem_parameter_key <- function(lhs, op, rhs) {
  paste(lhs, op, rhs)
}

.grsem_parameter_type <- function(pt) {
  out <- rep("other", nrow(pt))
  out[pt$op == "=~"] <- "loadings"
  out[pt$op == "~"] <- "regressions"
  out[pt$op == "~~" & pt$lhs == pt$rhs] <- "variances"
  out[pt$op == "~~" & pt$lhs != pt$rhs] <- "residual_covariances"
  out
}

.grsem_assign_groups <- function(pt_free, groups, penalize, ungrouped) {
  keys <- .grsem_parameter_key(pt_free$lhs, pt_free$op, pt_free$rhs)
  labels <- ifelse(nzchar(pt_free$label), pt_free$label, keys)
  eligible <- .grsem_parameter_type(pt_free) %in% penalize
  assigned <- rep(NA_character_, nrow(pt_free))

  if (is.null(groups)) {
    assigned[eligible] <- paste0("singleton:", labels[eligible])
  } else if (is.atomic(groups) && !is.null(names(groups))) {
    for (nm in names(groups)) {
      hit <- which(labels == nm | keys == nm)
      if (!length(hit)) {
        stop("group definition did not match parameter: ", nm, call. = FALSE)
      }
      if (any(!is.na(assigned[hit]))) {
        stop("a parameter was assigned to more than one group", call. = FALSE)
      }
      assigned[hit] <- as.character(groups[[nm]])
    }
    assigned[!eligible] <- NA_character_
  } else if (is.data.frame(groups) && "group" %in% names(groups)) {
    selectors <- intersect(c("label", "lhs", "op", "rhs"), names(groups))
    if (!length(selectors)) {
      stop("data-frame groups need at least one selector column", call. = FALSE)
    }
    for (row in seq_len(nrow(groups))) {
      hit <- rep(TRUE, nrow(pt_free))
      for (selector in selectors) {
        value <- groups[[selector]][row]
        if (!is.na(value) && nzchar(as.character(value))) {
          target <- if (selector == "label") labels else pt_free[[selector]]
          hit <- hit & target == as.character(value)
        }
      }
      idx <- which(hit & eligible)
      if (any(!is.na(assigned[idx]))) {
        stop("overlapping groups are not supported in version 0.1", call. = FALSE)
      }
      assigned[idx] <- as.character(groups$group[row])
    }
  } else {
    stop("groups must be NULL, a named vector, or a selector data frame",
         call. = FALSE)
  }

  l1_only <- eligible & is.na(assigned) & ungrouped == "l1"
  list(groups = assigned, l1_only = l1_only, labels = labels, keys = keys)
}

.grsem_build_spec <- function(backend, groups, penalize,
                              penalty_factor, group_weights, ungrouped) {
  fit <- backend$fit
  pt <- lavaan::parameterTable(fit)
  free_ids <- sort(unique(pt$free[pt$free > 0L]))
  pt_free <- pt[match(free_ids, pt$free), , drop = FALSE]
  full_spec <- .grsem_assign_groups(pt_free, groups, penalize, ungrouped)
  mapping <- list(map = backend$full_to_opt, scale = backend$full_scale,
                  n_opt = length(backend$x_unpenalized))

  opt_groups <- rep(NA_character_, mapping$n_opt)
  opt_l1_only <- rep(FALSE, mapping$n_opt)
  opt_labels <- rep(NA_character_, mapping$n_opt)
  for (i in seq_len(nrow(pt_free))) {
    j <- mapping$map[i]
    if (j == 0L) next
    proposed <- full_spec$groups[i]
    if (!is.na(proposed) && !is.na(opt_groups[j]) && opt_groups[j] != proposed) {
      stop("equality-constrained parameters must share one group", call. = FALSE)
    }
    if (!is.na(proposed)) opt_groups[j] <- proposed
    opt_l1_only[j] <- opt_l1_only[j] || full_spec$l1_only[i]
    if (is.na(opt_labels[j])) opt_labels[j] <- full_spec$labels[i]
  }
  for (j in seq_len(mapping$n_opt)) {
    members <- which(mapping$map == j)
    assigned <- !is.na(full_spec$groups[members]) | full_spec$l1_only[members]
    if (any(assigned) && !all(assigned)) {
      stop("all parameters in an equality class must share the same penalty status",
           call. = FALSE)
    }
  }

  default_pf <- as.numeric(!is.na(opt_groups) | opt_l1_only)
  if (is.null(penalty_factor)) {
    pf <- default_pf
  } else if (!is.null(names(penalty_factor))) {
    pf <- default_pf
    for (nm in names(penalty_factor)) {
      hit <- which(opt_labels == nm)
      if (!length(hit)) stop("penalty_factor name did not match: ", nm,
                             call. = FALSE)
      pf[hit] <- penalty_factor[[nm]]
    }
  } else if (length(penalty_factor) == mapping$n_opt) {
    pf <- as.numeric(penalty_factor)
  } else {
    stop("penalty_factor must be named or match the optimization dimension",
         call. = FALSE)
  }
  if (any(!is.finite(pf)) || any(pf < 0)) {
    stop("penalty_factor must be finite and non-negative", call. = FALSE)
  }
  grouped_idx <- which(!is.na(opt_groups) & pf > 0)
  if (any(!pf[grouped_idx] %in% c(1))) {
    stop("version 0.1 supports zero/one penalty factors within groups; use group_weights for group scaling",
         call. = FALSE)
  }

  penalized_full <- vapply(seq_len(nrow(pt_free)), function(i) {
    j <- mapping$map[i]
    j > 0 && pf[j] > 0 && !is.na(full_spec$groups[i])
  }, logical(1))
  sizes <- table(full_spec$groups[penalized_full])
  default_weights <- sqrt(as.numeric(sizes))
  names(default_weights) <- names(sizes)
  if (is.null(group_weights)) {
    gw <- default_weights
  } else {
    gw <- default_weights
    if (is.null(names(group_weights))) {
      if (length(group_weights) != length(gw)) {
        stop("unnamed group_weights has the wrong length", call. = FALSE)
      }
      gw[] <- group_weights
    } else {
      unknown <- setdiff(names(group_weights), names(gw))
      if (length(unknown)) stop("unknown group_weights: ", paste(unknown, collapse = ", "),
                                call. = FALSE)
      gw[names(group_weights)] <- group_weights
    }
  }
  if (any(!is.finite(gw)) || any(gw <= 0)) {
    stop("group_weights must be finite and positive", call. = FALSE)
  }

  l1_multiplier <- vapply(seq_len(mapping$n_opt), function(j) {
    members <- which(mapping$map == j)
    if (!length(members)) return(1)
    sum(abs(mapping$scale[members]))
  }, numeric(1))
  list(
    groups = opt_groups,
    l1_only = opt_l1_only,
    penalty_factor = pf,
    l1_weights = pf * l1_multiplier,
    group_weights = gw,
    labels = opt_labels,
    full_parameter_table = pt,
    pt_free = pt_free,
    full_to_opt = mapping$map
  )
}

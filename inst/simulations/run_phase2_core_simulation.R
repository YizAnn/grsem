source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
# Phase 2 pre-registered core simulation. Each core scenario uses 100
# replications. Extended scenarios are listed in phase2_simulation_design.csv and
# are not silently removed when computationally infeasible.
suppressPackageStartupMessages({
  library(grsem)
  library(lavaan)
})

B <- as.integer(Sys.getenv("GRSEM_PHASE2_REPS", "100"))
seed0 <- 20260722L
out <- output_dir("simulations/results_phase2")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

design <- data.frame(
  scenario = c("reflective_continuous", "reflective_binary_low_rate",
               "role_mimic", "role_problem_factor",
               "unequal_groups", "wrong_groups", "n100", "n500", "n1000",
               "p_near_n", "v1_v3_imbalance", "cross_version_homogeneous",
               "cross_version_heterogeneous", "partial_pooling_small_version"),
  tier = c(rep("core", 4), rep("extended", 10)),
  replications_planned = B,
  status = c(rep("scheduled", 4), rep("not_run_phase2_compute_boundary", 10)),
  stringsAsFactors = FALSE
)
write.csv(design, file.path(out, "phase2_simulation_design.csv"), row.names = FALSE)

candidate_names <- paste0("T", 1:6)
group_truth <- setNames(rep(paste0("G", 1:3), each = 2), candidate_names)
loading_truth <- setNames(c(0.65, 0.65, 0.55, 0, 0, 0), candidate_names)
group_definition <- setNames(unname(group_truth),
                             paste("PQ =~", candidate_names))

generate_reflective <- function(n, binary = FALSE, seed) {
  set.seed(seed)
  pq <- rnorm(n)
  dat <- data.frame(
    Q13 = 0.8 * pq + rnorm(n, sd = 0.6),
    Q15 = 0.75 * pq + rnorm(n, sd = 0.65)
  )
  thresholds <- setNames(seq(0.8, 1.5, length.out = 6), candidate_names)
  for (item in candidate_names) {
    latent <- loading_truth[item] * pq + rnorm(n)
    dat[[item]] <- if (binary) as.integer(latent > thresholds[item]) else latent
  }
  list(data = dat, thresholds = thresholds)
}

fit_recovery <- function(dat, binary, scenario, replicate, alpha,
                         threshold_truth = NULL) {
  start <- proc.time()[[3L]]
  method <- c("Group", "Sparse Group", "Lasso")[match(alpha, c(0, 0.5, 1))]
  ans <- try({
    fit <- grsem(
      paste0("PQ =~ Q13 + Q15 + ", paste(candidate_names, collapse = " + ")),
      dat, groups = group_definition, penalize = "loadings",
      rho = c(0.25, 0), alpha = alpha,
      backend = if (binary) "lavaan_dwls" else "lavaan_ml",
      ordered = if (binary) candidate_names else NULL,
      control = grsem_control(
        algorithm = "monotone_fista", max_iter = 1200,
        tolerance = 1e-4, approximate_tolerance = 5e-4, n_starts = 1
      )
    )
    row <- 1L
    pt <- summary(fit, row)$parameter_table
    est <- setNames(pt$est[pt$lhs == "PQ" & pt$op == "=~" &
                              pt$rhs %in% candidate_names],
                    pt$rhs[pt$lhs == "PQ" & pt$op == "=~" &
                             pt$rhs %in% candidate_names])
    est <- est[candidate_names]
    selected <- abs(est) > fit$control$zero_tolerance
    true_item <- loading_truth != 0
    selected_group <- tapply(selected, group_truth, any)
    true_group <- tapply(true_item, group_truth, any)
    threshold_rmse <- NA_real_
    if (binary) {
      pt0 <- summary(fit, 2L)$parameter_table
      threshold_est <- setNames(pt0$est[pt0$op == "|" &
                                           pt0$lhs %in% candidate_names],
                                pt0$lhs[pt0$op == "|" &
                                          pt0$lhs %in% candidate_names])
      threshold_est <- threshold_est[candidate_names]
      threshold_rmse <- sqrt(mean((threshold_est - threshold_truth)^2))
    }
    data.frame(
      scenario = scenario, replicate = replicate, method = method,
      backend = fit$backend$estimator, status = fit$path$status[row],
      lambda_max = fit$path$lambda_max[row], rho = fit$path$rho[row],
      lambda = fit$path$lambda[row],
      group_tpr = mean(selected_group[true_group]),
      group_fpr = mean(selected_group[!true_group]),
      parameter_tpr = mean(selected[true_item]),
      parameter_fpr = mean(selected[!true_item]),
      loading_rmse = sqrt(mean((est - loading_truth)^2)),
      threshold_rmse = threshold_rmse,
      nonzero_groups = fit$path$nonzero_groups[row],
      nonzero_parameters = fit$path$nonzero_penalized[row],
      prox_gradient_max = fit$path$prox_gradient_max[row],
      objective = fit$path$objective[row], iterations = fit$path$iterations[row],
      runtime_seconds = proc.time()[[3L]] - start, message = fit$path$message[row],
      stringsAsFactors = FALSE
    )
  }, silent = TRUE)
  if (inherits(ans, "try-error")) data.frame(
    scenario = scenario, replicate = replicate, method = method,
    backend = if (binary) "DWLS" else "ML", status = "failed",
    lambda_max = NA, rho = 0.25, lambda = NA, group_tpr = NA, group_fpr = NA,
    parameter_tpr = NA, parameter_fpr = NA, loading_rmse = NA,
    threshold_rmse = NA,
    nonzero_groups = NA, nonzero_parameters = NA, prox_gradient_max = NA,
    objective = NA, iterations = NA, runtime_seconds = proc.time()[[3L]] - start,
    message = as.character(ans), stringsAsFactors = FALSE
  ) else ans
}

role_models <- function(items) list(
  reflective = paste0("PQ =~ Q13 + Q15 + ", paste(items, collapse = " + "),
                      "\nSAT =~ S1 + S2 + S3\nSAT ~ PQ"),
  mimic = paste0("PQ =~ Q13 + Q15\nSAT =~ S1 + S2 + S3\nSAT ~ PQ\nPQ ~ ",
                 paste(paste0("M", 1:3), collapse = " + ")),
  problem_factor = paste0("PQ =~ Q13 + Q15\nSAT =~ S1 + S2 + S3\n",
                          "PROB =~ ", paste(items, collapse = " + "),
                          "\nPQ ~ PROB\nSAT ~ PQ + PROB")
)

generate_role <- function(n, role, seed) {
  set.seed(seed)
  module <- matrix(rnorm(n * 3), n, 3)
  problem <- rnorm(n)
  pq <- if (role == "mimic") -0.35 * module[, 1] - 0.25 * module[, 2] +
    rnorm(n, sd = 0.8) else if (role == "problem_factor")
      -0.45 * problem + rnorm(n, sd = 0.9) else rnorm(n)
  sat <- 0.65 * pq + if (role == "problem_factor") -0.25 * problem else 0 +
    rnorm(n, sd = 0.7)
  dat <- data.frame(
    Q13 = 0.8 * pq + rnorm(n, sd = 0.6), Q15 = 0.75 * pq + rnorm(n, sd = 0.65),
    S1 = 0.8 * sat + rnorm(n, sd = 0.6), S2 = 0.75 * sat + rnorm(n, sd = 0.65),
    S3 = 0.7 * sat + rnorm(n, sd = 0.7),
    M1 = module[, 1], M2 = module[, 2], M3 = module[, 3]
  )
  for (j in 1:6) {
    g <- ceiling(j / 2)
    source <- if (role == "problem_factor") problem else if (role == "mimic")
      module[, g] else pq
    dat[[paste0("T", j)]] <- 0.65 * source + rnorm(n, sd = 0.75)
  }
  dat
}

fit_role <- function(role, replicate) {
  dat <- generate_role(300, role, seed0 + 50000L + replicate +
                         if (role == "mimic") 1000L else 2000L)
  set.seed(seed0 + replicate)
  test <- sample(seq_len(nrow(dat)), 100)
  train <- dat[-test, ]; test_data <- dat[test, ]
  models <- role_models(candidate_names)
  losses <- setNames(rep(NA_real_, length(models)), names(models))
  statuses <- setNames(rep("failed", length(models)), names(models))
  for (name in names(models)) {
    fit <- try(sem(models[[name]], train, estimator = "ML", std.lv = TRUE,
                   fixed.x = FALSE, se = "none", test = "none"), silent = TRUE)
    if (!inherits(fit, "try-error") && lavInspect(fit, "converged")) {
      vars <- intersect(colnames(lavInspect(fit, "implied")$cov), names(test_data))
      sigma <- lavInspect(fit, "implied")$cov[vars, vars, drop = FALSE]
      losses[name] <- sqrt(mean((cov(test_data[, vars]) - sigma)^2))
      statuses[name] <- "converged"
    }
  }
  chosen <- if (all(!is.finite(losses))) NA_character_ else
    names(which.min(losses))
  data.frame(
    scenario = paste0("role_", role), replicate = replicate,
    true_role = if (role == "mimic") "mimic" else "problem_factor",
    selected_role = chosen, role_recovered = chosen ==
      if (role == "mimic") "mimic" else "problem_factor",
    reflective_loss = losses["reflective"], mimic_loss = losses["mimic"],
    problem_factor_loss = losses["problem_factor"],
    reflective_status = statuses["reflective"], mimic_status = statuses["mimic"],
    problem_factor_status = statuses["problem_factor"], stringsAsFactors = FALSE
  )
}

recovery_file <- file.path(out, "phase2_recovery_raw.csv")
role_file <- file.path(out, "phase2_role_raw.csv")
recovery_rows <- list(); role_rows <- list(); rid <- 0L
for (b in seq_len(B)) {
  for (binary in c(FALSE, TRUE)) {
    scenario <- if (binary) "reflective_binary_low_rate" else
      "reflective_continuous"
    generated <- generate_reflective(if (binary) 500 else 200, binary,
                                     seed0 + b + if (binary) 10000L else 0L)
    for (alpha in c(0, 0.5, 1)) {
      rid <- rid + 1L
      recovery_rows[[rid]] <- fit_recovery(
        generated$data, binary, scenario, b, alpha, generated$thresholds
      )
    }
  }
  role_rows[[2 * b - 1L]] <- fit_role("mimic", b)
  role_rows[[2 * b]] <- fit_role("problem_factor", b)
  write.csv(do.call(rbind, recovery_rows), recovery_file, row.names = FALSE)
  write.csv(do.call(rbind, role_rows), role_file, row.names = FALSE)
  if (b %% 5L == 0L) message("completed phase2 core replicate ", b, "/", B)
}

recovery <- do.call(rbind, recovery_rows)
roles <- do.call(rbind, role_rows)
summary <- do.call(rbind, lapply(split(recovery,
  interaction(recovery$scenario, recovery$method, drop = TRUE)), function(z) {
    ok <- z$status != "failed"
    data.frame(
      scenario = z$scenario[1], method = z$method[1], replications = nrow(z),
      failure_rate = mean(!ok), group_tpr = mean(z$group_tpr[ok], na.rm = TRUE),
      group_fpr = mean(z$group_fpr[ok], na.rm = TRUE),
      parameter_tpr = mean(z$parameter_tpr[ok], na.rm = TRUE),
      parameter_fpr = mean(z$parameter_fpr[ok], na.rm = TRUE),
      loading_rmse = mean(z$loading_rmse[ok], na.rm = TRUE),
      threshold_rmse = mean(z$threshold_rmse[ok], na.rm = TRUE),
      runtime_seconds = mean(z$runtime_seconds[ok], na.rm = TRUE)
    )
  }))
role_summary <- aggregate(role_recovered ~ scenario, roles,
                          function(x) mean(x, na.rm = TRUE))
write.csv(summary, file.path(out, "phase2_recovery_summary.csv"), row.names = FALSE)
write.csv(role_summary, file.path(out, "phase2_role_summary.csv"), row.names = FALSE)
writeLines(c(paste0("seed_base=", seed0), paste0("replications=", B),
             "core_scenarios=fixed before execution",
             "extended_scenarios=retained as not run, not deleted",
             "selection_threshold=rho 0.25; lambda_max recomputed per sample"),
           file.path(out, "phase2_reproducibility.txt"))
print(summary)
print(role_summary)

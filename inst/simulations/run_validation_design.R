source(file.path(Sys.getenv("GRSEM_PROJECT_ROOT"), "tools", "bootstrap.R"))
# Reproducible milestone-1 validation design for the grsem research prototype.
# Set GRSEM_SIM_REPS to increase replications (default 3).
suppressPackageStartupMessages({
  library(grsem)
  library(lavaan)
})

reps <- as.integer(Sys.getenv("GRSEM_SIM_REPS", "3"))
seed0 <- 20260718L
out_dir <- output_dir("simulations/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scenarios <- list(
  n200 = list(n=200, sizes=c(3,3,3), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0),
  n500 = list(n=500, sizes=c(3,3,3), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0),
  n1000 = list(n=1000, sizes=c(3,3,3), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0),
  sparse_within = list(n=500, sizes=c(4,4,4), signal=0.65, sparse=TRUE, wrong=FALSE, corr=0),
  high_indicator_correlation = list(n=500, sizes=c(3,3,3), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0.75),
  unequal_groups = list(n=500, sizes=c(2,5,8), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0),
  incorrect_groups = list(n=500, sizes=c(3,3,3), signal=0.65, sparse=TRUE, wrong=TRUE, corr=0),
  weak_loadings = list(n=500, sizes=c(3,3,3), signal=0.25, sparse=FALSE, wrong=FALSE, corr=0),
  strong_loadings = list(n=500, sizes=c(3,3,3), signal=0.85, sparse=FALSE, wrong=FALSE, corr=0),
  p_near_n = list(n=35, sizes=c(10,10,10), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0),
  p_greater_n = list(n=25, sizes=c(10,10,10), signal=0.65, sparse=FALSE, wrong=FALSE, corr=0)
)

make_case <- function(s, seed) {
  set.seed(seed)
  n <- s$n
  pq <- rnorm(n)
  sat <- 0.6 * pq + sqrt(1 - 0.6^2) * rnorm(n)
  dat <- data.frame(
    Q13 = 0.8 * pq + rnorm(n, sd=sqrt(1-0.8^2)),
    Q15 = 0.75 * pq + rnorm(n, sd=sqrt(1-0.75^2)),
    S1 = 0.8 * sat + rnorm(n, sd=sqrt(1-0.8^2)),
    S2 = 0.75 * sat + rnorm(n, sd=sqrt(1-0.75^2)),
    S3 = 0.7 * sat + rnorm(n, sd=sqrt(1-0.7^2))
  )
  candidate <- character()
  true_loading <- numeric()
  true_group <- character()
  module_shock <- if (s$corr > 0) {
    common <- rnorm(n)
    lapply(seq_along(s$sizes), function(g)
      s$corr * common + sqrt(1-s$corr^2) * rnorm(n))
  } else rep(list(rep(0, n)), length(s$sizes))
  for (g in seq_along(s$sizes)) {
    for (j in seq_len(s$sizes[g])) {
      nm <- paste0("T", g, "_", j)
      active_group <- g %in% c(1, 3)
      active_item <- active_group && (!s$sparse || j == 1L)
      a <- if (active_item) s$signal else 0
      residual_sd <- sqrt(max(0.05, 1-a^2-0.10*s$corr^2))
      dat[[nm]] <- a * pq + sqrt(0.10) * s$corr * module_shock[[g]] +
        rnorm(n, sd=residual_sd)
      candidate <- c(candidate, nm)
      true_loading[nm] <- a
      true_group[nm] <- paste0("G", g)
    }
  }
  assigned <- true_group
  if (s$wrong) assigned <- sample(assigned)
  model <- paste0(
    "PQ =~ Q13 + Q15 + ", paste(candidate, collapse=" + "), "\n",
    "SAT =~ S1 + S2 + S3\nPQ ~~ SAT"
  )
  groups <- setNames(unname(assigned), paste("PQ =~", names(assigned)))
  list(data=dat, model=model, groups=groups, truth=true_loading,
       true_group=true_group, assigned=assigned)
}

select_solution <- function(fit) {
  ok <- fit$path$converged & is.finite(fit$path$prox_gradient_max) &
    fit$path$prox_gradient_max <= 5e-4
  if (!any(ok)) return(NA_integer_)
  fit$path$solution[ok][which.min(fit$path$relative_bic_proxy[ok])]
}

one_method <- function(case, alpha, scenario, rep, weight_scheme="sqrt_size") {
  start <- proc.time()[[3]]
  ans <- tryCatch({
    weights <- if (weight_scheme == "sqrt_size") NULL else
      grsem_group_weights(case$model, case$data, case$groups,
                          method=weight_scheme)
    fit <- grsem(
      case$model, case$data, groups=case$groups, penalize="loadings",
      lambda=NULL, alpha=alpha, group_weights=weights,
      control=grsem_control(max_iter=1000, tolerance=1e-4,
        lambda_length=10, lambda_min_ratio=0.03, n_starts=2)
    )
    solution <- select_solution(fit)
    if (is.na(solution)) stop("no solution met the proximal-gradient threshold")
    row <- match(solution, fit$path$solution)
    pt <- summary(fit, solution)$parameter_table
    est <- setNames(pt$est[pt$lhs=="PQ" & pt$op=="=~" &
                            pt$rhs %in% names(case$truth)],
                    pt$rhs[pt$lhs=="PQ" & pt$op=="=~" &
                             pt$rhs %in% names(case$truth)])
    est <- est[names(case$truth)]
    selected <- abs(est) > fit$control$zero_tolerance
    true_active <- case$truth != 0
    group_selected <- tapply(selected, case$true_group, any)
    group_true <- tapply(true_active, case$true_group, any)
    refit <- try(post_selection_refit(fit, solution), silent=TRUE)
    fm <- if (inherits(refit, "try-error") || !lavInspect(refit,"converged"))
      rep(NA_real_,3) else fitMeasures(refit,c("cfi","rmsea","srmr"))
    corr_est <- if (inherits(refit, "try-error")) NA_real_ else {
      pe <- parameterEstimates(refit, standardized=TRUE)
      z <- pe[pe$lhs=="PQ" & pe$op=="~~" & pe$rhs=="SAT","std.all"]
      if (length(z)) z[1] else NA_real_
    }
    data.frame(
      scenario=scenario, replicate=rep,
      method=c("Group","Sparse-Group","Lasso")[match(alpha,c(0,.5,1))],
      alpha=alpha, weight_scheme=weight_scheme, n=nrow(case$data),
      p_candidates=length(case$truth), status="ok",
      group_tpr=mean(group_selected[group_true]),
      group_fpr=mean(group_selected[!group_true]),
      parameter_tpr=mean(selected[true_active]),
      parameter_fpr=mean(selected[!true_active]),
      loading_rmse=sqrt(mean((est-case$truth)^2)),
      loading_bias=mean(est-case$truth), latent_correlation=corr_est,
      latent_correlation_error=corr_est-0.6,
      cfi=unname(fm[1]), rmsea=unname(fm[2]), srmr=unname(fm[3]),
      raw_fit=fit$path$raw_fit[row], prox_gradient_max=fit$path$prox_gradient_max[row],
      iterations=fit$path$iterations[row], runtime_seconds=proc.time()[[3]]-start,
      selected_groups=fit$path$nonzero_groups[row],
      selected_parameters=fit$path$nonzero_penalized[row], message=""
    )
  }, error=function(e) data.frame(
    scenario=scenario, replicate=rep,
    method=c("Group","Sparse-Group","Lasso")[match(alpha,c(0,.5,1))],
    alpha=alpha, weight_scheme=weight_scheme, n=nrow(case$data),
    p_candidates=length(case$truth), status="failed",
    group_tpr=NA,group_fpr=NA,parameter_tpr=NA,parameter_fpr=NA,
    loading_rmse=NA,loading_bias=NA,latent_correlation=NA,
    latent_correlation_error=NA,cfi=NA,rmsea=NA,srmr=NA,raw_fit=NA,
    prox_gradient_max=NA,iterations=NA,runtime_seconds=proc.time()[[3]]-start,
    selected_groups=NA,selected_parameters=NA,message=conditionMessage(e)
  ))
  ans
}

rows <- list(); id <- 0L
for (scenario in names(scenarios)) for (rep in seq_len(reps)) {
  case <- make_case(scenarios[[scenario]], seed0 + rep + 100*match(scenario,names(scenarios)))
  for (alpha in c(0, 0.5, 1)) {
    id <- id+1L
    rows[[id]] <- one_method(case, alpha, scenario, rep)
  }
  message("completed ", scenario, " replicate ", rep)
}

# Direct scaling comparison in one stable scenario; this is separated from the
# recovery design so weighting and penalty type are not confounded.
case_scale <- make_case(scenarios$n500, seed0+9999L)
for (scheme in c("none","sqrt_size","hessian")) {
  id <- id+1L
  rows[[id]] <- one_method(case_scale, 0, "weight_scaling", 1L, scheme)
}

raw <- do.call(rbind, rows)
write.csv(raw, file.path(out_dir,"validation_raw.csv"), row.names=FALSE)
numeric_metrics <- c("group_tpr","group_fpr","parameter_tpr","parameter_fpr",
  "loading_rmse","loading_bias","latent_correlation_error","cfi","rmsea",
  "srmr","prox_gradient_max","runtime_seconds")
ok <- raw[raw$status=="ok",]
summary_rows <- if (nrow(ok)) do.call(rbind, lapply(
  split(ok, interaction(ok$scenario,ok$method,ok$weight_scheme,drop=TRUE)),
  function(z) {
    vals <- vapply(numeric_metrics, function(v) mean(z[[v]],na.rm=TRUE), numeric(1))
    data.frame(scenario=z$scenario[1],method=z$method[1],
      weight_scheme=z$weight_scheme[1],successful=nrow(z),t(vals),check.names=FALSE)
  })) else data.frame()
write.csv(summary_rows, file.path(out_dir,"validation_summary.csv"), row.names=FALSE)
write.csv(aggregate(status ~ scenario + method, raw, function(z) mean(z=="failed")),
  file.path(out_dir,"validation_failure_rates.csv"), row.names=FALSE)
writeLines(c(
  paste0("seed_base=",seed0), paste0("replications=",reps),
  "selection_rule=minimum relative BIC proxy among converged solutions with prox-gradient <= 5e-4",
  "ordinary likelihood-ratio inference was not assigned to penalized fits"
), file.path(out_dir,"reproducibility.txt"))

print(summary_rows)
print(aggregate(status ~ scenario + method, raw, function(z) mean(z=="failed")))

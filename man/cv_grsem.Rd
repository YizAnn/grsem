\name{cv_grsem}
\alias{cv_grsem}
\title{Sample-splitting evaluation}
\usage{cv_grsem(model, data, groups = NULL, penalize = "loadings",
rho = NULL, alpha = c(0, 0.5, 1), penalty_factor = NULL,
group_weights = NULL, estimator = "ML", std.lv = TRUE,
ungrouped = "unpenalized", meanstructure = FALSE,
control = grsem_control(), folds = 5L, seed = 1L, ...)}
\arguments{
\item{model,data,groups,penalize,rho,alpha,penalty_factor,group_weights,estimator,std.lv,ungrouped,meanstructure,control}{Arguments passed to \code{grsem}.}
\item{folds}{Number of folds (not a nested tuning pipeline).}
\item{seed}{Fold seed.}
\item{...}{Additional model arguments.}
}
\value{One row per fold and candidate, with test_nll, optimizer status, evaluation_status, failure_stage and message. Failed folds are retained with missing loss; no aggregation silently drops them.}
\description{Performs selection and unpenalized refitting inside each training split and evaluates multivariate-normal loss only on the held-out split.}

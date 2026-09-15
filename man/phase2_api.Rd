\name{phase2_api}
\alias{grsem_model}
\alias{fit_penalized}
\alias{bootstrap_penalized}
\alias{learn_candidate_roles}
\alias{grsem_prior}
\alias{fit_bayes}
\alias{grsem_bayes}
\alias{prior_predict}
\alias{posterior_predict}
\alias{compare_models}
\alias{print.grsem_model_spec}
\alias{summary.grsem_bayes_fit}
\title{Unified Phase 2 grsem research interfaces}
\description{
Defines one model specification shared by penalized estimation, auxiliary
candidate-role diagnostics, and explicitly experimental Bayesian backends.
The interfaces do not convert latent variables to factor scores for SEM
estimation.
}
\usage{
grsem_model(model, data, module_map = NULL, ordered = NULL,
  constraints = list(), metadata = list())

fit_penalized(spec, groups = NULL,
  backend = if (length(spec$ordered)) "lavaan_dwls" else "lavaan_ml", ...)

bootstrap_penalized(spec, groups, rho, alpha = c(0, 0.5, 1),
  B = 100L, seed = 1L, ...)

learn_candidate_roles(model, selected_items,
  module_map = model$module_map, constraints = model$constraints,
  methods = c("pc", "fci", "ges", "notears"), bootstrap = 100L,
  alpha = 0.05, m_max = 3L, seed = 1L)

grsem_prior(family = c("weak", "regularized_horseshoe",
  "group_horseshoe", "group_local"), global_scale = 0.1,
  slab_scale = 1, rope = c(-0.05, 0.05), sensitivity = list())

fit_bayes(spec, priors = grsem_prior(), backend = "cmdstanr",
  stan_model = c("continuous", "binary_probit", "hierarchical_versions"),
  chains = 4L, parallel_chains = chains, seed = 1L, ...)

grsem_bayes(spec, priors = grsem_prior(), backend = "cmdstanr", ...)

prior_predict(object, ...)

posterior_predict(object, variable = "y_rep", ...)

compare_models(fits, test_data = NULL, common_variables = NULL)

\method{print}{grsem_model_spec}(x, ...)

\method{summary}{grsem_bayes_fit}(object, ...)
}
\arguments{
\item{model, spec}{A lavaan syntax string or a \code{grsem_model_spec}.}
\item{data}{Analysis data frame.}
\item{module_map}{Named item-to-module map.}
\item{ordered}{Ordered/binary indicator names.}
\item{constraints}{Theory tiers, core-node definitions, and forbidden adjacencies.}
\item{metadata}{Provenance or backend-specific data.}
\item{groups}{Non-overlapping parameter groups.}
\item{rho}{Relative regularization path.}
\item{alpha}{Sparse-group mixing value or CI test level as documented by the function.}
\item{B, bootstrap}{Number of resamples.}
\item{seed}{Random seed.}
\item{selected_items}{Training-selected topic items.}
\item{methods}{Candidate structure-learning algorithms.}
\item{m_max}{Maximum conditioning-set size.}
\item{family, global_scale, slab_scale, rope, sensitivity}{Prior settings.}
\item{priors}{A \code{grsem_prior_spec}.}
\item{backend}{Estimation backend.}
\item{stan_model}{Experimental Stan program.}
\item{chains, parallel_chains}{MCMC chain counts.}
\item{object, x}{A fitted Bayesian object or model specification.}
\item{variable}{Posterior predictive generated quantity.}
\item{fits}{Named fitted-model list.}
\item{test_data}{Optional held-out observations.}
\item{common_variables}{Exact common variables required in every fitted model and test set.}
\item{...}{Additional backend arguments.}
}
\value{A model specification, fitted object, diagnostics, or comparison table.}
\details{
PC/FCI/GES, categorical DWLS and CmdStan interfaces are Experimental research prototypes. NOTEARS,
automatic general lavaan-to-Stan translation, and prior-predictive execution are
reported as not implemented rather than silently approximated.
}

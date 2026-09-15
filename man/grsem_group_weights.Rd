\name{grsem_group_weights}
\alias{grsem_group_weights}
\title{Construct group weights}
\usage{grsem_group_weights(model, data, groups, penalize = "loadings",
method = c("sqrt_size", "none", "hessian"), estimator = "ML",
std.lv = TRUE, meanstructure = FALSE, ...)}
\arguments{
\item{model,data,groups,penalize,estimator,std.lv,meanstructure}{Model arguments as in \code{grsem}.}
\item{method}{No scaling, square-root group size, or local Hessian scaling.}
\item{...}{Additional lavaan arguments.}
}
\value{Named positive group weights.}
\description{Computes alternative scaling schemes. Hessian weights use local curvature at the unpenalized solution and are normalized to geometric mean one; they are not effective degrees of freedom.}

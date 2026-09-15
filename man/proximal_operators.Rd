\name{proximal_operators}
\alias{soft_threshold}
\alias{group_soft_threshold}
\alias{sparse_group_prox}
\title{Proximal operators used by grsem}
\usage{
soft_threshold(x, threshold)

group_soft_threshold(x, threshold)

sparse_group_prox(x, groups, step_lambda, alpha, group_weights,
  penalty_factor, l1_only = rep(FALSE, length(x)))
}
\arguments{
\item{x}{Numeric vector.}
\item{threshold}{Nonnegative threshold.}
\item{groups}{Non-overlapping group ID per coordinate.}
\item{step_lambda}{Product of step size and lambda.}
\item{alpha}{Mixing value in zero to one.}
\item{group_weights}{Named group weights.}
\item{penalty_factor}{Per-coordinate L1 factors; zero is unchanged.}
\item{l1_only}{Ungrouped coordinates receiving only L1 shrinkage.}
}
\value{The thresholded numeric vector.}
\description{Closed-form soft, group-soft and non-overlapping sparse-group proximal maps.}

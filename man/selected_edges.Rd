\name{selected_edges}
\alias{selected_edges}
\title{Export selected structural regression edges}
\usage{selected_edges(object, solution, zero_tolerance = NULL)}
\arguments{
\item{object}{A \code{grsem_fit}.}
\item{solution}{One explicit solution identifier. Failed solutions are rejected.}
\item{zero_tolerance}{Finite nonnegative threshold; defaults to the fit control.
Only absolute exported estimates strictly greater than this are included.}
}
\description{Exports only lavaan \code{op == "~"} parameters as rhs-to-lhs edges.
Includes nonzero fixed and unpenalized regressions, identified in metadata.
Excludes measurement loadings, variances and covariances.}
\value{Data frame with from, to, estimate, parameter_id, lhs, op, rhs, label,
free, is_fixed, penalized, group, solution, status and zero_tolerance.}
\details{This is a graph of selected structural regressions, not causal discovery.
Thresholds use exported parameter units, including equality-constraint scaling.}

\name{selected_groups}
\alias{selected_groups}
\title{Selected parameter groups}
\usage{selected_groups(object, solution, tolerance = NULL)}
\arguments{
\item{object}{A \code{grsem_fit}.}
\item{solution}{Solution ID.}
\item{tolerance}{Optional zero threshold.}
}
\value{A data frame of group sizes, norms and selections.}
\description{Summarizes group selection at one regularization-path solution.}

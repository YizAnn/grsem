\name{post_selection_refit}
\alias{post_selection_refit}
\title{Post-selection unpenalized SEM refit}
\usage{post_selection_refit(object, solution, data = object$data,
store = FALSE, ...)}
\arguments{
\item{object}{A \code{grsem_fit}.}
\item{solution}{Solution ID.}
\item{data}{Training data.}
\item{store}{Return an updated object with the lavaan refit stored.}
\item{...}{Additional lavaan arguments.}
}
\value{A lavaan refit, or an updated \code{grsem_fit}.}
\description{Fixes zeroed parameters to zero and refits the selected structure without a penalty. Ordinary p-values remain conditional on the training-sample selection.}

\details{Only the Core ML backend is supported. Binary DWLS refits are rejected rather than silently converted to ML.}

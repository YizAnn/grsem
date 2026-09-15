\name{grsem_fit_methods}
\alias{print.grsem_fit}
\alias{summary.grsem_fit}
\alias{coef.grsem_fit}
\alias{plot.grsem_fit}
\alias{predict.grsem_fit}
\title{Methods for grouped regularized SEM fits}
\usage{
\method{print}{grsem_fit}(x, ...)
\method{summary}{grsem_fit}(object, solution = NULL, ...)
\method{coef}{grsem_fit}(object, solution = NULL, ...)
\method{plot}{grsem_fit}(x, type = c("path", "diagnostic"), ...)
\method{predict}{grsem_fit}(object, newdata = NULL, solution = NULL,
type = "lv", ...)
}
\arguments{
\item{x,object}{A \code{grsem_fit}.}
\item{solution}{Solution ID.}
\item{type}{Plot or prediction type.}
\item{newdata}{Optional new data for a stored post-selection refit.}
\item{...}{Additional method arguments.}
}
\value{Printed output, summaries, coefficients, a plot, or predictions.}
\description{Extraction and display methods. Prediction is deliberately available only after a post-selection lavaan refit has been stored.}

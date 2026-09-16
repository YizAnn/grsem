\name{grsem}
\alias{grsem}
\title{Grouped regularized structural equation modeling}
\usage{
grsem(model, data, groups = NULL, penalize = c("loadings"),
  lambda = NULL, rho = NULL, alpha = c(0, 0.25, 0.5), penalty_factor = NULL,
  group_weights = NULL, estimator = "ML", std.lv = TRUE,
  warm_start = TRUE, ungrouped = c("unpenalized", "l1"),
  meanstructure = FALSE, backend = c("lavaan_ml", "lavaan_dwls"),
  ordered = NULL, control = grsem_control(), ...)
}
\arguments{
\item{model}{lavaan model syntax.}
\item{data}{Complete numeric data frame.}
\item{groups}{Named parameter-to-group vector or selector data frame.}
\item{penalize}{Eligible parameter classes.}
\item{lambda}{Deprecated absolute path. Prefer \code{rho}.}
\item{rho}{Relative path \eqn{\lambda/\lambda_{max}}; recomputed for each sample.}
\item{alpha}{Mixing values: zero is Group Lasso and one is Lasso.}
\item{penalty_factor}{Per-coordinate factors; zero is unpenalized.}
\item{group_weights}{Positive named group weights.}
\item{estimator}{Only \code{"ML"} in version 0.1.}
\item{std.lv}{Passed to lavaan.}
\item{warm_start}{Use the preceding path solution as a start.}
\item{ungrouped}{Treatment of eligible but ungrouped parameters.}
\item{meanstructure}{Passed to lavaan.}
\item{backend}{Validated continuous ML or experimental ordered DWLS backend.}
\item{ordered}{Binary indicators for the experimental DWLS backend.}
\item{control}{Optimizer controls.}
\item{...}{Additional arguments for \code{lavaan::sem}.}
}
\value{A \code{grsem_fit} research-prototype object.}
\description{Fits non-overlapping Group Lasso and Sparse Group Lasso paths on a lavaan-exported continuous-variable, complete-data ML objective. All free SEM parameters are jointly estimated; only the selected coordinates receive proximal shrinkage.}
\examples{
# An overidentified CFA with four informative indicators and independent errors.
set.seed(1)
n <- 400
f <- rnorm(n)
d <- data.frame(
  x1 = 0.8 * f + rnorm(n, sd = 0.8),
  x2 = 0.9 * f + rnorm(n, sd = 0.8),
  x3 = 0.7 * f + rnorm(n, sd = 0.8),
  x4 = 0.85 * f + rnorm(n, sd = 0.8)
)
g <- c("f =~ x2"="candidate", "f =~ x3"="candidate",
       "f =~ x4"="candidate")
fit <- grsem("f =~ x1 + x2 + x3 + x4", d, groups=g,
             rho=c(0.2, 0), alpha=0,
             control=grsem_control(max_iter=100, tolerance=1e-3,
                                   n_starts=1))
fit
}

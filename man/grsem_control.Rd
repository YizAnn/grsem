\name{grsem_control}
\alias{grsem_control}
\title{Optimizer controls}
\usage{grsem_control(algorithm = c("monotone_fista", "fista", "pg"),
max_iter = 2000L, tolerance = 1e-5, approximate_tolerance = 1e-4,
objective_tolerance = 1e-8, initial_step = 1, backtrack = 0.5,
min_step = 1e-12, max_backtrack = 60L, lambda_length = 30L,
lambda_min_ratio = 0.01, zero_tolerance = 1e-6, n_starts = 2L,
restart = TRUE, checkpoint_path = NULL, checkpoint_every = 100L,
resume = FALSE, keep_history = FALSE)}
\arguments{
\item{algorithm}{Proximal gradient, FISTA, or monotone/restarted FISTA.}
\item{max_iter}{Maximum proximal-gradient iterations.}
\item{tolerance}{Maximum gradient-mapping tolerance.}
\item{approximate_tolerance}{First-order tolerance for approximately converged solutions.}
\item{objective_tolerance}{Relative objective-change tolerance.}
\item{initial_step}{Initial step size.}
\item{backtrack}{Step multiplier after a failed majorization check.}
\item{min_step}{Smallest step.}
\item{max_backtrack}{Maximum line-search attempts.}
\item{lambda_length}{Automatic path length.}
\item{lambda_min_ratio}{Automatic path endpoint relative to lambda max.}
\item{zero_tolerance}{Selection threshold.}
\item{n_starts}{Deterministic starts per solution.}
\item{restart}{Restart FISTA when momentum is counterproductive.}
\item{checkpoint_path}{Optional checkpoint directory.}
\item{checkpoint_every}{Checkpoint interval.}
\item{resume}{Resume compatible checkpoints.}
\item{keep_history}{Retain iteration histories.}
}
\value{A control list.}
\description{Creates controls for the backtracking proximal-gradient solver.}

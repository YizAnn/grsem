# Method and implementation map

## Research question and workflow

Can a user-specified SEM be regularized at the parameter-group and individual
parameter levels while preserving explicit separation between selection,
structural graph interpretation, refitting, and evaluation?

SEM specification → structured regularization → parameter selection →
structural edge extraction → DAG validity diagnostics → post-selection refit →
evaluation.

The code is a research prototype built around lavaan, not a new causal discovery
algorithm. All distributed executable examples use synthetic observations.

## Objective and parameter mapping

The Core objective is the lavaan-exported smooth loss plus:

~~~text
lambda * [(1-alpha) * sum_g w_g ||x_g||_2 + alpha * sum_j v_j |x_j|]
~~~

R/groups.R maps the free parameter table into optimization coordinates. The
default eligible class is loadings; regressions, variances and residual
covariances can be selected explicitly. Only regressions can become directed
structural edges. Fixed parameters are not optimization coordinates.

Groups are non-overlapping. Eligible ungrouped coordinates are unpenalized by
default, or L1-only if requested. Zero penalty factors bypass shrinkage while
remaining part of the smooth optimization. Supported simple equality classes
share coordinates, group membership and penalty status.

Default group weights use the square root of the number of penalized full
parameter-table entries in each group. This need not equal the number of
independent optimization coordinates under equality constraints. L1 weights
also incorporate equality-map scales.

## Paths and optimization

- lambda controls total regularization.
- alpha=0 gives group lasso; alpha=1 gives L1; intermediate values give
  sparse-group lasso.
- rho=lambda/lambda_max is a relative setting, not an absolute penalty.
- R/optimizer.R computes lambda_max from the gradient at a restricted fit
  with penalized coordinates set to zero. This is not a global-optimality
  certificate for a nonconvex SEM.
- Each alpha and each training/bootstrap sample recomputes lambda_max.
- R/prox.R applies coordinate soft thresholding followed by group thresholding.
- R/optimizer.R combines smooth-gradient updates, backtracking, box bounds,
  warm starts and optional deterministic alternative starts.
- PG, FISTA and monotone-FISTA branches exist. The public example uses PG.
  The existing restart expression may clear momentum frequently; no
  acceleration-rate or speedup guarantee is made.

Diagnostics distinguish converged, approximately_converged and failed.
Convergence uses a proximal-gradient mapping and objective change, rather
than merely the fact that an iteration loop terminated.

## Parameter selection and structural edges

R/grsem.R retains path diagnostics and parameter coordinates. selected_groups()
in R/methods.R marks a group selected if any relevant coordinate exceeds the
zero threshold; this is not a group-norm significance test.

selected_edges() in R/graph.R takes one explicit nonfailed solution and exports
only op "~" as rhs → lhs. Loadings (=~), variances and covariances (~~) never
become directed edges through this interface.

The output includes from, to, estimate, parameter_id, lhs, op, rhs, label, free,
is_fixed, penalized, group, solution, status and zero_tolerance. Nonzero fixed
and unpenalized regressions are included and identified in metadata.

Only abs(estimate) > zero_tolerance is retained, using exported parameter units.
At equality with the threshold, the parameter is excluded. The default is
the fitted control's zero_tolerance. Approximate fits retain their status;
failed fits cannot be exported.

## DAG legality

validate_dag() uses Kahn topological sorting and iterative DFS for a cycle
witness. Duplicate edges are deduplicated; a self-loop is a cycle. An empty
edge set is a DAG. Additional isolated nodes can be supplied explicitly.

The result has is_dag, topological_order, cycle and nodes. A non-DAG has no
topological order; its cycle repeats the starting node at the end.

This checks a selected regression graph. It does not identify causes,
orient unidentified edges, learn a new SEM, or implement NOTEARS. The SEM
numerical-validity checks are separate from graph acyclicity.

## Post-selection refit

R/refit.R fixes near-zero penalized coordinates in the original parameter
table and calls lavaan ML without the penalty. It does not translate an
arbitrary graph into SEM syntax. Experimental DWLS refits are rejected.

The refit uses the original optimizer-coordinate zero rule; selected_edges()
uses exported parameter units. Equality-map scaling can make near-threshold
decisions differ. Standard refit p-values do not adjust for the selection step.

## Evaluation and failure accounting

R/cv.R fits and refits inside each training fold and evaluates Gaussian NLL on
the held-out fold. It aligns named covariance variables, means and test columns;
missing model means become named zero vectors.

Every fold and alpha/rho candidate remains in the result. Optimizer status is
separate from evaluation_status. Failed fits, refits and NLL/chol computations
have an explicit failure_stage, message and missing loss. Downstream code
must not silently discard these failures when comparing candidates.

R/compare_models.R retains successful and failed models in a fixed schema.
It reports fit indices and optional covariance reconstruction error on exactly
the requested common variables. It does not mechanically compare ordinary
information criteria across different observed variable sets.

CV is not a complete nested model-selection procedure. Relative AIC/BIC/EBIC
fields are proxies. This project does not claim calibrated selective inference.

## Implementation map

| Step | Files and main functions |
|---|---|
| Specification | R/model_spec.R: grsem_model |
| Fit wrapper / path | R/penalized_fit.R: fit_penalized; R/grsem.R: grsem |
| Backend / groups | R/backend_lavaan.R; R/groups.R; R/weights.R |
| Proximal optimization | R/prox.R; R/optimizer.R |
| Selection / edge export | R/methods.R: selected_groups; R/graph.R: selected_edges |
| DAG diagnosis | R/graph.R: validate_dag |
| Refit | R/refit.R: post_selection_refit |
| Evaluation | R/cv.R: cv_grsem; R/compare_models.R: compare_models |
| Synthetic executable example | inst/examples/quick_start.R |

Bayesian and PC/FCI/GES are separate Experimental interfaces. NOTEARS and
prior_predict are Not Implemented. See [Limitations](LIMITATIONS.md).

# Limitations and capability boundaries

## Core scope

The Core workflow is a continuous, complete-data ML research prototype for a
user-specified SEM. Non-overlapping parameter groups and a restricted class of
simple equality constraints are supported. This is not a general production
SEM platform; missing-data/FIML support, broad categorical inference and all
possible constraints are not established.

The objective may be nonconvex. Converged or approximately_converged is a
numerical status, not a proof of global optimality, identification, or a
scientifically correct model. Existing FISTA restart behavior can frequently
remove momentum, so acceleration is not guaranteed.

## Selection, graph validity, and inference

Threshold-based sparsity is not a significance test. Small changes in scaling
or zero_tolerance can change near-zero edges. Only structural regressions
enter selected_edges; measurement loadings and covariances do not.

validate_dag checks graph legality, not causal validity. It cannot establish
the absence of confounding, orient a partially directed graph, discover a
causal structure, or prove statistical identification.

Refitting is ML-only and uses the original parameter table, not a generic
graph-to-SEM translation. Refit zero rules use optimization coordinates,
while exported edge thresholds use original parameter units; equality-map
scales matter near the threshold. Ordinary refit p-values are not
selection-adjusted inference.

## Evaluation

cv_grsem retains failed candidates and folds. Do not silently remove failed
folds or average away missing losses. Strict convergence, approximate
convergence, and evaluation success are different states.

The workflow does not automatically perform nested tuning followed by an
independent final evaluation. Information-criterion fields on regularized
fits remain relative proxies. Comparing models using different observed
variables requires care; requested covariance-comparison variables must
exist in every fitted model and held-out set.

The synthetic example uses a fixed demonstration rho, not an optimized
setting. Its noise edge is not always eliminated. It is a reproducibility
example, not a perfect-recovery, superiority, or business-impact claim.

## Experimental and unimplemented interfaces

| Interface | Status | Important boundary |
|---|---|---|
| Binary DWLS | Experimental | Limited existing backend; no general ordered-data guarantee |
| Bayesian / CmdStan | Experimental | Explicit Stan data required; no general lavaan-to-Stan conversion |
| posterior_predict | Experimental | Extracts generated draws; not arbitrary new-data prediction |
| PC / FCI / GES | Experimental | External-package candidate-role wrappers |
| NOTEARS | Not Implemented | Returns a placeholder status even if external software exists |
| prior_predict | Not Implemented | Stops; no dedicated prior predictive program |

Bayesian prior settings are recorded but are not fully wired or validated
against Stan data. Divergence-only status reporting is insufficient for
sampling convergence. Its benchmark entry is disabled in the Core workflow;
no Bayesian result is claimed here.

Candidate-role wrappers use module averages rather than fitted latent scores.
Forbidden adjacencies affect skeleton search; theory tiers are used for
post-hoc direction annotation, not general search-direction constraints.
Matrix-to-direction interpretation and statistical coverage remain
experimental. These wrappers are not connected automatically to the Core
selected-edge/refit pipeline.

## Environment and output policy

Recorded local checks use Windows 11, R 4.4.3 and lavaan 0.7-2.
Linux/macOS, the minimum declared R 4.3, clean-machine dependency installation,
and online CI have not been verified. No PDF manual build was performed.

R >= 4.3 and lavaan >= 0.7-2 are requirements; the availability of compatible
packages from a chosen repository must be confirmed at installation time.
The dependency setup script rejects an older lavaan.

Core repository scripts load the candidate installation from .Rlib and write
under the configured project root, with path checks. This is an application
policy, not an operating-system sandbox for arbitrary R code. Unchanged
Experimental Bayesian internals are outside the Core output guarantee.

Use the staging helper for R CMD build when temp files are inside the
repository. Dependencies, logs, generated artifacts, staging and historical
construction documents are not public source inputs.

## Data and attribution

All distributed executable examples generate synthetic data. Restricted
survey observations and results derived from them are not distributed.
The project does not make real-world performance claims from synthetic runs.

Zhuoer Zhang originated the research direction and method ideas, with no
collaborators or academic supervisors. The project uses established methods
and third-party software; dependency authorship and algorithmic novelty are
not claimed. The repository retains an MIT license declaration.

# grsem

**An R research prototype for structured regularization of user-specified structural equation models, with explicit selection, graph diagnostics, and evaluation.**

**SEM specification → structured regularization → parameter selection → structural edge extraction → DAG validity diagnostics → post-selection refit → evaluation**

The Core workflow connects statistical modeling, proximal optimization, graph
validation, and failure-aware evaluation. A valid DAG here is a property of
selected structural regressions—not evidence of causal discovery.

## Research motivation

The project was initially inspired by customer satisfaction SEM research:
when candidate parameters have meaningful groups, can selection preserve that
structure while keeping estimation, graph interpretation, and evaluation
separate?

grsem explores this question with group and sparse-group penalties over SEM
parameters. This repository contains **synthetic examples only**; no restricted
survey data or results derived from such data are distributed.

## Core workflow

1. Specify candidate relationships in lavaan syntax with grsem_model().
2. Fit a structured regularization path with fit_penalized().
3. Choose an explicit solution and inspect its numerical status.
4. Extract nonzero structural regressions with selected_edges().
5. Diagnose directed cycles with validate_dag().
6. Fix selected zeros and refit the original parameter table without a penalty.
7. Evaluate with fold-level Gaussian NLL or comparable held-out covariance error.

Selection does not invent new candidate edges. Measurement loadings and
covariances are not exported as directed structural edges.

## Method overview

For backend loss F and free optimization coordinates x, the objective is:

~~~text
F(x) + lambda * [(1 - alpha) * sum_g w_g ||x_g||_2
                + alpha * sum_j v_j |x_j|]
~~~

- alpha=0: group lasso; alpha=1: L1; intermediate values: sparse-group lasso.
- rho=lambda/lambda_max defines a relative path; lambda_max is recomputed for
  each alpha and training/bootstrap sample.
- Parameters with zero penalty factors remain estimated but bypass shrinkage.
- PG and FISTA variants use proximal updates and backtracking. The example
  uses PG; acceleration limitations are documented rather than hidden.
- selected_edges() includes only op "~", directed rhs → lhs. It includes
  nonzero fixed/unpenalized regressions with metadata, and applies a strict
  absolute-estimate threshold.
- validate_dag() returns is_dag, a topological order for a DAG, or a concrete
  directed cycle otherwise. It is not NOTEARS or a structure-learning method.

See [Method](docs/METHOD.md) for the implementation map and
[Limitations](docs/LIMITATIONS.md) for interpretation boundaries.

## Quick Start

Requires **R >= 4.3** and **lavaan >= 0.7-2** with lav_export_estimation().
A lower lavaan version is not a compatible substitute. This project is not
claimed to be a CRAN release.

From the repository root, put a supported R installation on PATH and verify
Rscript.exe --version. Example PowerShell setup uses no machine-specific paths:

~~~powershell
$env:GRSEM_PROJECT_ROOT = (Get-Location).Path
$env:R_LIBS_USER = Join-Path $env:GRSEM_PROJECT_ROOT '.Rlib'
$env:TMP = Join-Path $env:GRSEM_PROJECT_ROOT 'artifacts/tmp'
$env:TEMP = $env:TMP
New-Item -ItemType Directory -Force -Path $env:R_LIBS_USER,$env:TMP | Out-Null
Rscript.exe --vanilla tools/setup_core.R
R.exe CMD INSTALL --library=.Rlib .
Rscript.exe --vanilla -e "source('tools/bootstrap.R'); source('inst/examples/quick_start.R')"
~~~

Alternatively, configure R_LIBS to reuse a compatible dependency library.
Always install this candidate grsem into its own .Rlib; the bootstrap verifies
the loaded package location. A clean-machine dependency installation has not
yet been validated; package availability must satisfy the version requirements.

The complete runnable example is [quick_start.R](inst/examples/quick_start.R).
Its essential modeling steps are:

~~~r
library(grsem)
set.seed(20260915)
n <- 400L
x <- rnorm(n)
z <- rnorm(n)
m <- 0.8 * x + rnorm(n, sd = 0.7)
y <- 0.7 * m + 0.25 * x + rnorm(n, sd = 0.7)
dat <- data.frame(x = x, z = z, m = m, y = y)

spec <- grsem_model("m ~ x + z\ny ~ m + x + z", dat)
groups <- c("m ~ x" = "m_inputs", "m ~ z" = "m_inputs",
            "y ~ m" = "y_inputs", "y ~ x" = "y_inputs",
            "y ~ z" = "y_inputs")
fit <- fit_penalized(
  spec, groups = groups, penalize = "regressions",
  rho = c(0.2, 0), alpha = 0.5,
  control = grsem_control(algorithm = "pg", max_iter = 3000,
                          tolerance = 1e-5, n_starts = 1)
)
solution <- fit$path$solution[fit$path$rho == 0.2]
edges <- selected_edges(fit, solution, zero_tolerance = 1e-6)
dag <- validate_dag(edges, nodes = names(dat))
stopifnot(dag$is_dag)
refit <- post_selection_refit(fit, solution, data = dat)
~~~

rho=0.2 is a predeclared demonstration setting, **not a CV-selected optimum**.
The script checks optimizer status and refit convergence. Evaluation is a
separate step: cv_grsem() repeats selection/refitting within training folds
and returns held-out losses with explicit failure records.

## Example output

Recorded from a fresh synthetic run on Windows 11, R 4.4.3, lavaan 0.7-2:

| from | to | estimate |
|---|---|---:|
| x | m | 0.715051 |
| z | m | 0.012083 |
| m | y | 0.532072 |
| x | y | 0.282887 |

The rho=0.2 fit converged, is_dag was TRUE, the returned topological order was
x, z, m, y, and the ML refit converged. The noise path z→m remains nonzero:
this demonstrates a working pipeline, **not perfect graph recovery**.
Small numerical differences across environments are possible.

## Validation / Tests

Local validation on **2026-09-15**, Windows 11 / R 4.4.3 / lavaan 0.7-2:

- 38 test blocks, 110 passing assertions; no failures, warnings, or skips.
- Independent R --vanilla synthetic Quick Start completed.
- Source-package build, including the synthetic vignette, completed.
- R CMD check --no-manual: **Status: OK**, including package tests and vignette
  rebuilding. Optional Suggests were not forced.
- No speedup, production-readiness, or broad benchmark claim is made.

Tests cover numerical gradient agreement, proximal operators, parameter groups,
relative paths, failure preservation, Gaussian loss alignment, graph cycles,
edge semantics, and the end-to-end synthetic workflow.

~~~powershell
Rscript.exe --vanilla tools/run_core_tests.R
Rscript.exe --vanilla tools/stage_core.R
$sourcePath = (Get-Content artifacts/build-source-path.txt -Raw).Trim()
R.exe CMD build --no-manual $sourcePath
$env:_R_CHECK_FORCE_SUGGESTS_ = 'false'
New-Item -ItemType Directory -Force -Path artifacts/check | Out-Null
R.exe CMD check --no-manual --output=artifacts/check grsem_0.2.0.9000.tar.gz
~~~

Use staging when the R temporary directory is inside this repository; building
the root directly can recursively copy temporary files. Logs, installed
libraries, and generated results stay local and are excluded from publication.
The GitHub Actions configuration has not been executed on a remote repository.

## Capability status

| Capability | Status | Scope |
|---|---|---|
| Continuous, complete-data regularized SEM | Core | lavaan ML; user-specified model |
| Group / sparse-group penalty and relative paths | Core | Non-overlapping parameter groups |
| selected_edges / validate_dag | Core | Structural regressions and graph legality |
| ML post-selection refit and evaluation | Core | Failure records; no selective-inference guarantee |
| Binary DWLS | Experimental | Existing limited backend and regression test |
| Bayesian / posterior_predict | Experimental | CmdStan templates and draw extraction; explicit Stan data |
| PC / FCI / GES | Experimental | Auxiliary candidate-role wrappers, not the Core pipeline |
| NOTEARS | Not Implemented | Placeholder branch |
| prior_predict | Not Implemented | Reserved interface that stops explicitly |

## Limitations

SEM objectives can be nonconvex: numerical convergence does not establish a
global optimum. Graph validity does not establish causal identification.
Ordinary refit p-values do not account for selection, IC columns are relative
proxies, and CV is not an automatic nested tuning/final-evaluation pipeline.

Validation is local to the documented environment; Linux/macOS, R 4.3, clean
dependency installation, and online CI remain unverified. Read
[all limitations](docs/LIMITATIONS.md) before interpreting model results.

## Contribution statement

**Zhuoer Zhang** originated the research direction and methodological ideas.
There were **no collaborators or academic supervisors** on this project.

The implementation builds on lavaan and the R ecosystem; experimental wrappers
delegate to their respective external packages. This statement does not claim
authorship of those dependencies or novelty for established regularization and
graph algorithms.

Contact: [zhangze2048@ruc.edu.cn](mailto:zhangze2048@ruc.edu.cn).
Package author and maintainer: Zhuoer Zhang. License: MIT.

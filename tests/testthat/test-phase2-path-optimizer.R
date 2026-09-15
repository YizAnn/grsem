test_that("rho is remapped through a sample-specific lambda_max", {
  dat <- simulate_one_factor(400, seed = 31)
  scaled <- dat * 10
  rho <- c(0.5, 0.3, 0.2, 0.1)
  control <- grsem_control(max_iter = 1500, tolerance = 2e-5,
                           approximate_tolerance = 2e-4, n_starts = 1)
  fit1 <- grsem(candidate_model, dat, groups = candidate_groups,
                rho = rho, alpha = 0.5, control = control)
  fit2 <- grsem(candidate_model, scaled, groups = candidate_groups,
                rho = rho, alpha = 0.5, control = control)
  absolute <- grsem(candidate_model, scaled, groups = candidate_groups,
                    lambda = fit1$path$lambda, alpha = 0.5,
                    control = control)
  expect_equal(fit1$path$lambda, rho * fit1$path$lambda_max,
               tolerance = 1e-12)
  expect_equal(fit2$path$lambda, rho * fit2$path$lambda_max,
               tolerance = 1e-12)
  expect_equal(fit1$path$nonzero_penalized,
               fit2$path$nonzero_penalized)
  expect_true(any(absolute$path$nonzero_penalized !=
                    fit2$path$nonzero_penalized))
})

test_that("PG and accelerated solvers reach comparable objectives", {
  dat <- simulate_one_factor(300, seed = 44)
  algorithms <- c("pg", "fista", "monotone_fista")
  fits <- lapply(algorithms, function(algorithm) {
    grsem(candidate_model, dat, groups = candidate_groups,
          rho = 0.2, alpha = 0.5,
          control = grsem_control(
            algorithm = algorithm, max_iter = 2500, tolerance = 2e-5,
            approximate_tolerance = 2e-4, n_starts = 1
          ))
  })
  objectives <- vapply(fits, function(x) x$path$objective, numeric(1))
  expect_lt(max(objectives) - min(objectives), 2e-4)
  expect_true(all(vapply(fits, function(x)
    x$path$status %in% c("converged", "approximately_converged"), logical(1))))
})

test_that("status uses the Phase 2 three-level vocabulary", {
  dat <- simulate_one_factor(250, seed = 18)
  fit <- grsem(candidate_model, dat, groups = candidate_groups,
               rho = c(1, 0.2, 0), alpha = 0.5,
               control = grsem_control(max_iter = 500, tolerance = 1e-4,
                                       approximate_tolerance = 1e-3,
                                       n_starts = 1))
  expect_true(all(fit$path$status %in%
                    c("converged", "approximately_converged", "failed")))
})

test_that("experimental binary DWLS lambda zero matches lavaan", {
  set.seed(901)
  n <- 600
  eta <- rnorm(n)
  dat <- data.frame(
    y1 = as.integer(eta + rnorm(n) > 0.2),
    y2 = as.integer(0.8 * eta + rnorm(n) > 0.5),
    y3 = as.integer(0.7 * eta + rnorm(n) > 0.8),
    y4 = as.integer(rnorm(n) > 1.2)
  )
  model <- "f =~ y1 + y2 + y3 + y4"
  groups <- c("f =~ y2" = "signal", "f =~ y3" = "signal",
              "f =~ y4" = "noise")
  fit <- grsem(model, dat, groups = groups, rho = 0, alpha = 0.5,
               backend = "lavaan_dwls", ordered = names(dat),
               control = grsem_control(n_starts = 1))
  expect_equal(as.numeric(coef(fit, 1)),
               as.numeric(fit$backend$x_unpenalized), tolerance = 1e-8)
  expect_equal(fit$backend$estimator, "DWLS")
})

test_that("unified model specification is retained", {
  dat <- simulate_one_factor(200, seed = 71)
  spec <- grsem_model(candidate_model, dat,
                      module_map = setNames(c("a", "a"), c("x4", "x5")))
  fit <- fit_penalized(spec, groups = candidate_groups,
                       rho = 0.2, alpha = 0.5,
                       control = grsem_control(max_iter = 1000,
                                               tolerance = 1e-4,
                                               n_starts = 1))
  expect_s3_class(spec, "grsem_model_spec")
  expect_identical(fit$model_spec$schema_version, "0.2")
})

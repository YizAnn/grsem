test_that("Gaussian NLL aligns covariance, mean, and reordered test columns", {
  sigma <- matrix(c(2, .4, .4, 1), 2,
                  dimnames = list(c("x", "y"), c("x", "y")))
  data <- data.frame(y = c(-1, 2, .5), extra = 99, x = c(.3, 1, -2))
  mu <- c(y = .2, x = -.3)
  z <- sweep(as.matrix(data[c("x", "y")]), 2, mu[c("x", "y")], "-")
  reference <- mean(.5 * (2 * log(2*pi) + log(det(sigma)) +
                            rowSums((z %*% solve(sigma)) * z)))
  expect_equal(grsem:::.grsem_gaussian_nll(sigma, mu, data), reference)
  expect_equal(grsem:::.grsem_gaussian_nll(sigma, NULL, data),
               grsem:::.grsem_gaussian_nll(sigma, c(x = 0, y = 0), data))
  expect_equal(grsem:::.grsem_gaussian_nll(sigma, c(-.3, .2), data), reference)
  expect_error(grsem:::.grsem_gaussian_nll(sigma, c(x = 1), data), "mean")
  expect_error(grsem:::.grsem_gaussian_nll(sigma, NULL, data["x"]), "all model")
  bad <- sigma; bad[1, 1] <- -1
  expect_error(grsem:::.grsem_gaussian_nll(bad, NULL, data), "not positive")
})
test_that("CV records every candidate when a fold fit throws", {
  local_mocked_bindings(grsem = function(...) stop("deliberate fit failure"),
                        .package = "grsem")
  result <- cv_grsem("unused", data.frame(x = 1:12), folds = 3,
                     rho = c(.5, 0), alpha = c(0, 1))
  expect_equal(nrow(result), 12L)
  expect_true(all(result$evaluation_status == "failed"))
  expect_true(all(result$failure_stage == "fit"))
  expect_true(all(is.na(result$test_nll)))
})
test_that("CV records a real chol error and continues other folds", {
  dat <- simulate_one_factor(160, seed = 400)
  local_mocked_bindings(.grsem_test_nll = function(fit, newdata) {
    chol(matrix(c(1, 2, 2, 1), 2))
  }, .package = "grsem")
  result <- cv_grsem(candidate_model, dat, groups = candidate_groups,
                     rho = 0, alpha = 1, folds = 2)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$evaluation_status == "failed"))
  expect_true(all(result$failure_stage == "nll"))
  expect_true(all(grepl("not positive", result$message)))
})
test_that("CV catches refit errors separately", {
  local_mocked_bindings(post_selection_refit = function(...) stop("refit error"),
                        .package = "grsem")
  result <- cv_grsem(candidate_model, simulate_one_factor(120),
                     groups = candidate_groups, rho = 0, alpha = 1, folds = 2)
  expect_true(all(result$failure_stage == "refit"))
})
test_that("CV without meanstructure yields finite held-out scores", {
  result <- cv_grsem(candidate_model, simulate_one_factor(160, seed = 402),
                     groups = candidate_groups, rho = 0, alpha = 1, folds = 2)
  expect_equal(nrow(result), 2L)
  expect_true(all(result$evaluation_status == "ok"))
  expect_true(all(is.finite(result$test_nll)))
})
test_that("mixed valid and invalid model comparisons preserve all rows", {
  fit <- lavaan::cfa(candidate_model, simulate_one_factor(200), std.lv = TRUE)
  result <- compare_models(list(good = fit, bad = NULL))
  expect_identical(result$model, c("good", "bad"))
  expect_identical(result$status, c("converged", "failed"))
  expect_true(is.na(result$CFI[2]))
  expect_match(result$message[2], "not a fitted")
  expect_equal(nrow(compare_models(list())), 0L)
  bad_common <- compare_models(list(good = fit), test_data = simulate_one_factor(50),
                              common_variables = c("x1", "absent"))
  expect_identical(bad_common$status, "failed")
})
test_that("checkpoint outputs cannot escape the configured candidate root", {
  root <- Sys.getenv("GRSEM_PROJECT_ROOT")
  skip_if(!nzchar(root), "requires explicit candidate root")
  expect_error(grsem:::.grsem_output_path("../outside.rds"), "traversal")
  expect_error(grsem:::.grsem_output_path(file.path(dirname(root), "outside.rds")),
               "inside")
  expect_match(grsem:::.grsem_output_path("artifacts/safe.rds"), "artifacts/safe.rds",
               fixed = TRUE)
})

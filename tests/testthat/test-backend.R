test_that("lavaan export gradient agrees with numerical differentiation", {
  skip_if_not_installed("numDeriv")
  dat <- simulate_one_factor(300)
  backend <- grsem:::.grsem_prepare_backend(
    candidate_model, dat, "ML", TRUE, FALSE
  )
  x <- backend$x_unpenalized + seq_along(backend$x_unpenalized) * 1e-4
  numerical <- numDeriv::grad(backend$fn, x)
  analytic <- backend$gr(x)
  expect_equal(analytic, numerical, tolerance = 2e-5)
})

test_that("lambda zero is the unpenalized lavaan solution", {
  dat <- simulate_one_factor(300)
  fit <- grsem(
    candidate_model, dat, groups = candidate_groups,
    lambda = 0, alpha = c(0, 0.5, 1),
    control = grsem_control(n_starts = 1)
  )
  target <- fit$backend$x_unpenalized
  for (i in seq_len(nrow(fit$path))) {
    expect_equal(as.numeric(coef(fit, fit$path$solution[i])),
                 as.numeric(target), tolerance = 1e-8)
  }
})

test_that("unpenalized coordinates are not directly zeroed", {
  dat <- simulate_one_factor(300)
  fit <- grsem(
    candidate_model, dat, groups = candidate_groups,
    lambda = 0.2, alpha = 0.5,
    control = grsem_control(max_iter = 500, tolerance = 1e-4, n_starts = 1)
  )
  unpenalized <- fit$spec$penalty_factor == 0
  expect_true(all(abs(fit$coefficients[[1]][unpenalized]) > 0))
})

test_that("alpha one is invariant to the irrelevant group partition", {
  dat <- simulate_one_factor(300, seed = 7)
  g1 <- candidate_groups
  g2 <- setNames(paste0("singleton_", seq_along(candidate_groups)),
                 names(candidate_groups))
  control <- grsem_control(max_iter = 1000, tolerance = 1e-5, n_starts = 1)
  fit1 <- grsem(candidate_model, dat, groups = g1,
                lambda = 0.03, alpha = 1, control = control)
  fit2 <- grsem(candidate_model, dat, groups = g2,
                lambda = 0.03, alpha = 1, control = control)
  expect_equal(fit1$coefficients[[1]], fit2$coefficients[[1]], tolerance = 2e-4)
})

test_that("post-selection refit keeps penalized zeros fixed at zero", {
  dat <- simulate_one_factor(300, seed = 21)
  fit <- grsem(candidate_model, dat, groups = candidate_groups,
               lambda = 10, alpha = 1,
               control = grsem_control(max_iter = 1000, tolerance = 1e-4,
                                       n_starts = 1))
  refit <- post_selection_refit(fit, 1, se = "none", test = "none")
  pt <- lavaan::parameterTable(refit)
  zeroed <- pt$lhs == "f" & pt$op == "=~" &
    pt$rhs %in% c("x2", "x3", "x4", "x5") & pt$free == 0
  expect_true(any(zeroed))
  expect_equal(pt$est[zeroed], rep(0, sum(zeroed)))
})

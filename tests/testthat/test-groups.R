test_that("selector data frames define nonoverlapping groups", {
  dat <- simulate_one_factor(200)
  selectors <- data.frame(
    lhs = "f", op = "=~", rhs = c("x2", "x3", "x4", "x5"),
    group = c("signal", "signal", "noise", "noise")
  )
  fit <- grsem(candidate_model, dat, groups = selectors,
               lambda = 0, alpha = 0)
  expect_equal(sort(names(fit$spec$group_weights)), c("noise", "signal"))
})

test_that("overlapping selectors are rejected", {
  dat <- simulate_one_factor(200)
  selectors <- data.frame(
    lhs = c("f", "f"), op = c("=~", "=~"),
    rhs = c("x2", "x2"), group = c("a", "b")
  )
  expect_error(grsem(candidate_model, dat, groups = selectors,
                     lambda = 0, alpha = 0), "overlapping")
})

test_that("simple equality constraints are preserved", {
  dat <- simulate_one_factor(300, seed = 9)
  model <- "f =~ x1 + equal*x2 + equal*x3 + x4 + x5"
  fit <- grsem(model, dat, groups = c(equal = "signal"),
               lambda = c(0.02, 0), alpha = 0,
               control = grsem_control(max_iter = 800, tolerance = 1e-4,
                                       n_starts = 1))
  pt <- summary(fit, 1)$parameter_table
  equal_rows <- pt$label == "equal"
  expect_equal(pt$est[equal_rows][1], pt$est[equal_rows][2], tolerance = 1e-10)
})

test_that("fixed loadings remain fixed", {
  dat <- simulate_one_factor(300, seed = 11)
  model <- "f =~ 1*x1 + x2 + x3 + x4 + x5"
  fit <- grsem(model, dat, groups = candidate_groups,
               lambda = 0.03, alpha = 0, std.lv = FALSE,
               control = grsem_control(max_iter = 800, tolerance = 1e-4,
                                       n_starts = 1))
  pt <- summary(fit, 1)$parameter_table
  expect_equal(pt$est[pt$lhs=="f" & pt$op=="=~" & pt$rhs=="x1"], 1)
})

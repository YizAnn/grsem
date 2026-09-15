test_that("soft threshold has the expected closed form", {
  expect_equal(soft_threshold(c(-2, -0.5, 0, 0.5, 2), 1),
               c(-1, 0, 0, 0, 1))
})

test_that("group soft threshold has the expected closed form", {
  expect_equal(group_soft_threshold(c(3, 4), 2), c(1.8, 2.4))
  expect_equal(group_soft_threshold(c(3, 4), 5), c(0, 0))
})

test_that("sparse group prox is soft then group threshold", {
  x <- c(3, 4, 7)
  got <- sparse_group_prox(
    x, groups = c("a", "a", NA), step_lambda = 1, alpha = 0.5,
    group_weights = c(a = 2), penalty_factor = c(1, 1, 0)
  )
  expected_group <- group_soft_threshold(soft_threshold(c(3, 4), 0.5), 1)
  expect_equal(got, c(expected_group, 7))
})

test_that("singleton group lasso equals weighted lasso", {
  x <- c(-2, 0.5, 3)
  groups <- c("a", "b", "c")
  got <- sparse_group_prox(
    x, groups, step_lambda = 0.4, alpha = 0,
    group_weights = c(a = 1, b = 1, c = 1),
    penalty_factor = rep(1, 3)
  )
  expect_equal(got, soft_threshold(x, 0.4))
})

test_that("zero penalty factor is untouched by the proximal map", {
  got <- sparse_group_prox(
    c(1, 1), groups = c("a", "a"), step_lambda = 100, alpha = 0.5,
    group_weights = c(a = sqrt(1)), penalty_factor = c(1, 0)
  )
  expect_equal(got[2], 1)
})

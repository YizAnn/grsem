test_that("weight schemes are positive and named", {
  dat <- simulate_one_factor(250)
  w_size <- grsem_group_weights(candidate_model, dat, candidate_groups,
                                method = "sqrt_size")
  w_none <- grsem_group_weights(candidate_model, dat, candidate_groups,
                                method = "none")
  expect_equal(w_size, c(noise = sqrt(2), signal = sqrt(2)))
  expect_equal(w_none, c(noise = 1, signal = 1))
  skip_if_not_installed("numDeriv")
  w_h <- grsem_group_weights(candidate_model, dat, candidate_groups,
                             method = "hessian")
  expect_true(all(is.finite(w_h) & w_h > 0))
  expect_equal(exp(mean(log(w_h))), 1, tolerance = 1e-8)
})

edge_frame <- function(from = character(), to = character()) {
  data.frame(from = from, to = to, stringsAsFactors = FALSE)
}
expect_topological <- function(edges, result) {
  expect_true(result$is_dag)
  positions <- match(edges$from, result$topological_order)
  expect_true(all(positions < match(edges$to, result$topological_order)))
  expect_length(result$cycle, 0)
}
expect_cycle_edges <- function(edges, result) {
  expect_false(result$is_dag)
  expect_length(result$topological_order, 0)
  expect_identical(result$cycle[1], tail(result$cycle, 1))
  pairs <- paste(edges$from, edges$to)
  expect_true(all(paste(head(result$cycle, -1), tail(result$cycle, -1)) %in% pairs))
}
test_that("simple DAG has a valid topological order", {
  e <- edge_frame(c("a", "a", "b"), c("b", "c", "c"))
  expect_topological(e, validate_dag(e))
})
test_that("two-node directed cycle is returned explicitly", {
  e <- edge_frame(c("a", "b"), c("b", "a"))
  expect_cycle_edges(e, validate_dag(e))
})
test_that("three-node directed cycle is returned explicitly", {
  e <- edge_frame(c("a", "b", "c"), c("b", "c", "a"))
  expect_cycle_edges(e, validate_dag(e))
})
test_that("disconnected DAG and isolated vertices are retained", {
  e <- edge_frame(c("a", "c"), c("b", "d"))
  result <- validate_dag(e, nodes = "isolated")
  expect_topological(e, result)
  expect_setequal(result$topological_order, c("a", "b", "c", "d", "isolated"))
})
test_that("empty edge sets are DAGs", {
  result <- validate_dag(edge_frame())
  expect_true(result$is_dag)
  expect_length(result$topological_order, 0)
  expect_length(result$cycle, 0)
  expect_identical(validate_dag(edge_frame(), nodes = "a")$topological_order, "a")
})
test_that("duplicates, self-loops and invalid endpoints are handled", {
  expect_topological(edge_frame("a", "b"),
                     validate_dag(edge_frame(c("a", "a"), c("b", "b"))))
  expect_cycle_edges(edge_frame("a", "a"), validate_dag(edge_frame("a", "a")))
  expect_error(validate_dag(edge_frame("", "a")), "node names")
  expect_error(validate_dag(edge_frame(NA_character_, "a")), "node names")
  expect_error(validate_dag(data.frame(x = "a")), "from and to")
})

# A mixed parameter-table fixture isolates graph semantics from optimizer noise.
mixed_fit <- function() {
  pt <- data.frame(
    lhs = c("y", "f", "x", "y", "z", "w"),
    op = c("~", "=~", "~~", "~~", "~", "~"),
    rhs = c("x", "item", "x", "x", "y", "z"),
    label = c("beta", "", "", "", "fixed", "unpenalized"),
    free = c(1L, 2L, 3L, 4L, 0L, 5L), ustart = c(NA, NA, NA, NA, .5, NA)
  )
  structure(list(
    path = data.frame(solution = 7L, status = "converged"),
    coefficients = list(c(.2, .8, 1, .4, .3)),
    backend = list(get_coef = function(x) x),
    spec = list(full_parameter_table = pt, full_to_opt = 1:5,
                penalty_factor = c(1, 1, 0, 0, 0),
                groups = c("reg", "loading", NA, NA, NA)),
    control = grsem_control(zero_tolerance = .1)
  ), class = "grsem_fit")
}
test_that("only structural regressions export, with parameter metadata", {
  e <- selected_edges(mixed_fit(), 7L)
  expect_identical(e$from, c("x", "y", "z"))
  expect_identical(e$to, c("y", "z", "w"))
  expect_true(all(e$op == "~"))
  expect_identical(e$is_fixed, c(FALSE, TRUE, FALSE))
  expect_identical(e$penalized, c(TRUE, FALSE, FALSE))
  expect_equal(e$estimate, c(.2, .5, .3))
  expect_true(all(e$solution == 7L))
})
test_that("measurement loading does not enter the DAG", {
  e <- selected_edges(mixed_fit(), 7L)
  expect_false(any(e$from == "f" | e$to == "item"))
})
test_that("variance and residual covariance do not enter the DAG", {
  e <- selected_edges(mixed_fit(), 7L)
  expect_false(any(e$parameter_id %in% c(3L, 4L)))
})
test_that("zero threshold is strict and caller-controlled", {
  fit <- mixed_fit()
  fit$coefficients[[1]][1] <- -.2
  expect_false(any(selected_edges(fit, 7L, .2)$parameter_id == 1L))
  expect_true(any(selected_edges(fit, 7L, .199)$estimate == -.2))
  fit$coefficients[[1]][1] <- 0
  expect_false(any(selected_edges(fit, 7L, 0)$parameter_id == 1L))
  expect_equal(nrow(selected_edges(fit, 7L, 1)), 0L)
  expect_true(validate_dag(selected_edges(fit, 7L, 1))$is_dag)
  expect_error(selected_edges(fit, 7L, -1), "zero_tolerance")
  expect_error(selected_edges(fit, 99L), "unknown solution")
  fit$path$status <- "failed"
  expect_error(selected_edges(fit, 7L), "failed solution")
})

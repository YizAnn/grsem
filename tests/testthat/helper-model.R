simulate_one_factor <- function(n = 500L, seed = 1L) {
  set.seed(seed)
  eta <- rnorm(n)
  loading <- c(0.9, 0.8, 0.7, 0, 0)
  x <- vapply(loading, function(a) a * eta + rnorm(n, sd = 0.6), numeric(n))
  x <- as.data.frame(x)
  names(x) <- paste0("x", seq_along(loading))
  x
}

candidate_groups <- c(
  "f =~ x2" = "signal",
  "f =~ x3" = "signal",
  "f =~ x4" = "noise",
  "f =~ x5" = "noise"
)

candidate_model <- "f =~ x1 + x2 + x3 + x4 + x5"

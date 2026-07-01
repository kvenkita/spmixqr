test_that("resolve_weights normalizes to mean 1 and preserves raw sum", {
  rw <- resolve_weights(c(1, 2, 3, 4), "frequency", data = NULL, n = 4L)
  expect_equal(mean(rw$w), 1)
  expect_equal(rw$sum_raw, 10)
  expect_equal(rw$raw, c(1, 2, 3, 4))
  expect_true(rw$weighted)
  expect_identical(rw$type, "frequency")
})

test_that("resolve_weights NULL is inert (unit weights, not flagged weighted)", {
  rw <- resolve_weights(NULL, "sampling", data = NULL, n = 3L)
  expect_equal(rw$w, rep(1, 3))
  expect_null(rw$raw)
  expect_false(rw$weighted)
  expect_equal(rw$sum_raw, 3)
})

test_that("resolve_weights resolves a column name and a one-sided formula", {
  d <- data.frame(y = 1:3, wt = c(2, 2, 2))
  expect_equal(resolve_weights("wt", "precision", d, 3L)$raw, c(2, 2, 2))
  expect_equal(resolve_weights(~wt, "precision", d, 3L)$raw, c(2, 2, 2))
})

test_that("resolve_weights rejects bad input", {
  expect_error(resolve_weights(c(1, -1, 1), "sampling", NULL, 3L), "non-negative")
  expect_error(resolve_weights(c(1, 1), "sampling", NULL, 3L), "length")
  expect_error(resolve_weights(c(0, 0, 0), "sampling", NULL, 3L), "all")
  expect_error(resolve_weights(c(1, NA, 1), "sampling", NULL, 3L), "finite")
  expect_error(resolve_weights(~wt, "sampling", data = NULL, n = 3L), "data")
})

test_that("weighted gate matches row-duplication (frequency semantics)", {
  set.seed(11)
  n <- 60; Z <- cbind(1, rnorm(n)); G <- 2
  P <- spmixqr:::normalize_rows(matrix(runif(n * G), n, G))
  Pen <- diag(c(0, 1e-3))
  wct <- sample(1:3, n, replace = TRUE)
  fit_w <- spmixqr:::pen_irls_multinom(Z, P, Pen, w = wct)
  dup <- rep(seq_len(n), wct)
  fit_d <- spmixqr:::pen_irls_multinom(Z[dup, ], P[dup, ], Pen)
  expect_equal(as.numeric(fit_w$gamma), as.numeric(fit_d$gamma), tolerance = 1e-5)
})

test_that("unit weights reproduce the unweighted gate exactly", {
  set.seed(12)
  n <- 40; Z <- cbind(1, rnorm(n)); P <- spmixqr:::normalize_rows(matrix(runif(n * 2), n, 2))
  Pen <- diag(c(0, 1e-3))
  a <- spmixqr:::pen_irls_multinom(Z, P, Pen)
  b <- spmixqr:::pen_irls_multinom(Z, P, Pen, w = rep(1, n))
  expect_identical(a$gamma, b$gamma)
})

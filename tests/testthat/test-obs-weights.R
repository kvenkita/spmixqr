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

test_that("spatial_em_fit unit weights reproduce the unweighted fit", {
  set.seed(21)
  d <- sim_spmixqr(n = 120, G = 2, tau = 0.5, seed = 21)
  f0 <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
                variance = "none", control = spmixqr_control(nstart = 2L, seed = 1))
  f1 <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
                weights = rep(1, nrow(d$data)),
                variance = "none", control = spmixqr_control(nstart = 2L, seed = 1))
  expect_equal(f0$beta_const, f1$beta_const, tolerance = 1e-8)
  expect_equal(f0$loglik, f1$loglik, tolerance = 1e-8)
})

test_that("frequency weights match row duplication in the unpenalized regime", {
  ## Under mean-1 weight normalization, frequency = duplication holds EXACTLY only
  ## without a penalty (scale-invariant weighted quantile regression). G=1 with the
  ## spatial terms off isolates that regime.
  set.seed(31)
  d <- sim_spmixqr(n = 90, G = 2, tau = 0.5, seed = 31)
  wct <- sample(1:3, nrow(d$data), replace = TRUE)
  fw <- spmixqr(y ~ x, data = d$data, G = 1, tau = 0.5,
                spatial_coef = FALSE, spatial_gate = FALSE,
                weights = wct, weights_type = "frequency", variance = "none",
                control = spmixqr_control(nstart = 1L, seed = 7))
  dup <- rep(seq_len(nrow(d$data)), wct)
  fd <- spmixqr(y ~ x, data = d$data[dup, ], G = 1, tau = 0.5,
                spatial_coef = FALSE, spatial_gate = FALSE, variance = "none",
                control = spmixqr_control(nstart = 1L, seed = 7))
  expect_equal(fw$beta_const, fd$beta_const, tolerance = 1e-6)
})

test_that("fit stores weight metadata; unweighted stores NULL prior_weights", {
  d <- sim_spmixqr(n = 50, G = 2, tau = 0.5, seed = 5)
  f <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
               weights = runif(50, 0.5, 2), weights_type = "precision",
               variance = "none", control = spmixqr_control(nstart = 1L))
  expect_identical(f$weights_type, "precision")
  expect_equal(mean(f$weights), 1)
  expect_length(f$prior_weights, 50)
  f0 <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
                variance = "none", control = spmixqr_control(nstart = 1L))
  expect_null(f0$prior_weights)
})

test_that("sampling meat squares weights; frequency meat is linear; unit weights inert", {
  set.seed(41)
  n <- 80; Xt <- cbind(1, rnorm(n)); y <- Xt %*% c(1, 2) + rnorm(n)
  w <- runif(n, 0.3, 1); beta <- c(1, 2); Pen <- diag(c(0, 1e-3)); h <- 0.3
  ow <- runif(n, 0.5, 2)
  base <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h)
  same <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h, ow = rep(1, n), wtype = "sampling")
  expect_equal(base, same, tolerance = 1e-10)              # ow = 1 is a no-op
  samp <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h, ow = ow, wtype = "sampling")
  freq <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h, ow = ow, wtype = "frequency")
  expect_false(isTRUE(all.equal(samp, freq)))               # meat scaling differs
})

test_that("coef sandwich unit-weight inertness holds for the frequency branch too", {
  set.seed(43)
  n <- 60; Xt <- cbind(1, rnorm(n)); y <- Xt %*% c(1, 2) + rnorm(n)
  w <- runif(n, 0.3, 1); beta <- c(1, 2); Pen <- diag(c(0, 1e-3)); h <- 0.3
  base <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h)
  freq1 <- spmixqr:::coef_sandwich_vcov(Xt, y, 0.5, w, beta, Pen, h, ow = rep(1, n), wtype = "frequency")
  expect_equal(base, freq1, tolerance = 1e-10)
})

test_that("end-to-end sandwich runs with live weights for each type", {
  set.seed(44)
  d <- sim_spmixqr(n = 70, G = 2, tau = 0.5, seed = 44)
  wt <- runif(70, 0.5, 2)
  for (ty in c("sampling", "frequency", "precision")) {
    f <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
                 weights = wt, weights_type = ty, variance = "sandwich",
                 control = spmixqr_control(nstart = 1L))
    expect_false(is.null(f$vcov))
    expect_true(all(is.finite(diag(f$vcov$coef[[1]]))))
  }
})

test_that("weighted bootstrap runs and returns SEs for each type", {
  skip_on_cran()
  set.seed(61)
  d <- sim_spmixqr(n = 80, G = 2, tau = 0.5, seed = 61)
  wt <- runif(80, 0.5, 2)
  for (ty in c("sampling", "frequency", "precision")) {
    f <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
                 weights = wt, weights_type = ty, variance = "none",
                 control = spmixqr_control(nstart = 1L, boot_B = 20L, seed = 3))
    V <- spmixqr:::bootstrap_vcov(f, d$data)
    expect_true(is.list(V$coef) && length(V$coef) == 2)
    expect_true(all(is.finite(diag(V$coef[[1]]))))
  }
})

test_that("frequency bootstrap carries weights into each refit", {
  set.seed(62)
  d <- sim_spmixqr(n = 70, G = 2, tau = 0.5, seed = 62)
  wt <- sample(1:3, 70, replace = TRUE)
  f <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
               weights = wt, weights_type = "frequency", variance = "none",
               control = spmixqr_control(nstart = 1L, boot_B = 15L, seed = 9))
  expect_silent(spmixqr:::bootstrap_vcov(f, d$data))
})

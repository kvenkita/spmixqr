test_that("gate penalised gradient matches numDeriv and Newton matches optim", {
  skip_if_not_installed("numDeriv")
  set.seed(7); n <- 150
  Z <- cbind(1, matrix(rnorm(n * 2), n, 2)); q1 <- ncol(Z); G <- 3
  P <- spmixqr:::normalize_rows(matrix(runif(n * G), n, G))
  Pen <- diag(c(1e-3, 0.5, 0.5))
  g0 <- matrix(rnorm(q1 * (G - 1)) * 0.2, q1, G - 1)
  num <- numDeriv::grad(function(v)
    spmixqr:::gate_objective(matrix(v, q1, G - 1), Z, P, Pen), as.numeric(g0))
  pis <- spmixqr:::softmax_rows(Z %*% g0); ana <- numeric(q1 * (G - 1))
  for (a in seq_len(G - 1)) {
    r <- P[, a + 1] - pis[, a + 1]
    ana[((a - 1) * q1 + 1):(a * q1)] <- crossprod(Z, r) - Pen %*% g0[, a]
  }
  expect_lt(max(abs(num - ana)), 1e-5)
  fit <- spmixqr:::pen_irls_multinom(Z, P, Pen)
  opt <- optim(as.numeric(g0),
               function(v) -spmixqr:::gate_objective(matrix(v, q1, G - 1), Z, P, Pen),
               method = "BFGS")
  expect_equal(spmixqr:::gate_objective(fit$gamma, Z, P, Pen), -opt$value, tolerance = 1e-3)
})

test_that("smoothed penalised QR matches numDeriv and reduces to weighted_rq", {
  skip_if_not_installed("numDeriv")
  set.seed(7); n <- 250
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  y <- as.numeric(X %*% c(1, 2, -1) + (rexp(n) - rexp(n)) * 0.5)
  w <- runif(n, .3, 1); tau <- 0.5; h <- 0.2; Pen <- diag(c(0, 0, 0))
  b0 <- c(.5, 1, -.5)
  num <- numDeriv::grad(function(b)
    spmixqr:::coef_objective(b, X, y, tau, w, Pen, h), b0)
  e <- as.numeric(y - X %*% b0)
  ana <- -(crossprod(X, w * spmixqr:::psi_smooth(e, tau, h)) - Pen %*% b0)
  expect_lt(max(abs(num - as.numeric(ana))), 1e-5)
  sm <- spmixqr:::pen_smooth_wqr(X, y, tau, w, Pen, h = 0.02)
  wr <- mixqr::weighted_rq(X, y, tau, w = w)
  expect_lt(max(abs(sm$beta - wr)), 0.05)
})

test_that("roughness penalty shrinks spatial-slope coefficients", {
  set.seed(3); n <- 200
  X <- cbind(1, rnorm(n)); y <- as.numeric(X %*% c(0, 1) + rnorm(n) * 0.5)
  PenBig <- diag(c(0, 1e6))
  sm <- spmixqr:::pen_smooth_wqr(X, y, 0.5, rep(1, n), PenBig, h = 0.2)
  expect_lt(abs(sm$beta[2]), 1e-2)
})

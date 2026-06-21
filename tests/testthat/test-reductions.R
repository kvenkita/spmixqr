test_that("flat gate + flat slopes reduces to mixqr", {
  set.seed(42)
  d <- sim_spmixqr(n = 220, seed = 42)
  f0 <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                spatial_gate = FALSE, spatial_coef = FALSE, variance = "none",
                control = spmixqr_control(nstart = 3L, seed = 1))
  m <- mixqr::mixqr(y ~ x, data = d$data, m = 2, tau = 0.5)
  expect_equal(f0$loglik, as.numeric(stats::logLik(m)), tolerance = 1e-3)
})

test_that("G = 1 with no spatial terms matches quantreg::rq", {
  skip_if_not_installed("quantreg")
  set.seed(5); n <- 300
  dat <- data.frame(y = rnorm(n) + 2, x = rnorm(n))
  cc <- cbind(runif(n), runif(n))
  f1 <- spmixqr(y ~ x, dat, coords = cc, G = 1, tau = 0.5,
                spatial_gate = FALSE, spatial_coef = FALSE, variance = "none",
                control = spmixqr_control(nstart = 1L, seed = 1))
  rq <- quantreg::rq(y ~ x, tau = 0.5, data = dat)
  expect_equal(as.numeric(f1$beta_const), as.numeric(coef(rq)), tolerance = 1e-4)
})

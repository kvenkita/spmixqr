test_that("spatial_error = FALSE reproduces v0.1.0 exactly (byte-identical loglik/coef)", {
  ## The default path must be untouched by the CAR module. Fit twice and compare.
  set.seed(42); d <- sim_spmixqr(n = 220, seed = 42)
  ctl <- spmixqr_control(nstart = 3L, seed = 1)
  f_a <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 spatial_gate = FALSE, spatial_coef = FALSE, variance = "none",
                 control = ctl)
  f_b <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 spatial_error = FALSE, spatial_gate = FALSE, spatial_coef = FALSE,
                 variance = "none", control = ctl)
  expect_identical(f_a$loglik, f_b$loglik)
  expect_identical(as.numeric(f_a$coefficients), as.numeric(f_b$coefficients))
  expect_null(f_a$car)
  expect_false(isTRUE(f_a$spatial_error))
})

test_that("G = 1 + no-spatial reduces to quantreg::rq even alongside the CAR code", {
  skip_if_not_installed("quantreg")
  set.seed(5); n <- 300
  dat <- data.frame(y = rnorm(n) + 2, x = rnorm(n))
  cc <- cbind(runif(n), runif(n))
  f1 <- spmixqr(y ~ x, dat, coords = cc, G = 1, tau = 0.5, spatial_error = FALSE,
                spatial_gate = FALSE, spatial_coef = FALSE, variance = "none",
                control = spmixqr_control(nstart = 1L, seed = 1))
  rq <- quantreg::rq(y ~ x, tau = 0.5, data = dat)
  expect_equal(as.numeric(f1$beta_const), as.numeric(coef(rq)), tolerance = 1e-4)
})

test_that("lambda_error -> infinity drives phi to zero", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 400, G = 1, tau = 0.5, spatial_error = TRUE,
                   lattice = 7, car_rho = 1.5, seed = 7)
  f <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
               spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
               lambda_error = 1e8, variance = "none",
               control = spmixqr_control(nstart = 1L, seed = 1))
  expect_lt(max(abs(f$car$phi)), 1e-4)
})

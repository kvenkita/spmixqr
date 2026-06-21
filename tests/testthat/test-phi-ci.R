test_that("phi_surface ci/scale add SEs, credible flags, and exp() multiplier", {
  set.seed(1)
  d <- sim_spmixqr(n = 250, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 6, seed = 1)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "boot",
                 control = spmixqr_control(nstart = 1L, seed = 1, boot_B = 30L))
  ps <- phi_surface(fit, ci = TRUE, scale = "exp")
  expect_true(all(c("phi", "se", "lower", "upper", "credible", "mult") %in% names(ps)))
  expect_true(all(is.finite(ps$se)))
  expect_true(all(ps$se >= 0))
  expect_true(is.logical(ps$credible))
  expect_equal(ps$mult, exp(ps$phi), tolerance = 1e-10)
  expect_true(all(ps$lower <= ps$phi & ps$phi <= ps$upper))
  ## credible == interval excludes zero
  expect_equal(ps$credible, is.finite(ps$se) & (ps$lower > 0 | ps$upper < 0))
})

test_that("phi_surface(ci = TRUE) warns and returns NA SEs without a covariance", {
  set.seed(2)
  d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 5, seed = 2)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_warning(ps <- phi_surface(fit, ci = TRUE), "standard errors")
  expect_true(all(is.na(ps$se)))
})

test_that("phi standard errors equal the Tmat-propagated coefficient covariance", {
  set.seed(3)
  d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 5, seed = 3)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "sandwich", control = spmixqr_control(nstart = 1L, seed = 1))
  cb <- fit$car$car_block
  Tm <- as.matrix(fit$car$Tmat)
  Vred <- as.matrix(fit$vcov$coef[[1]])[cb, cb, drop = FALSE]
  se_manual <- sqrt(pmax(diag(Tm %*% Vred %*% t(Tm)), 0))
  expect_equal(phi_surface(fit, ci = TRUE)$se, se_manual, tolerance = 1e-8)
})

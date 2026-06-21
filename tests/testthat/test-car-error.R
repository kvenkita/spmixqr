test_that("G = 1 CAR spatial-error recovers beta and the phi surface", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 500, G = 1, tau = 0.5, spatial_error = TRUE,
                   lattice = 7, car_rho = 1.5, seed = 7)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_true(isTRUE(fit$spatial_error))
  ## phi recovery
  expect_gt(cor(fit$car$phi[, 1], d$truth$phi[, 1]), 0.85)
  ## beta recovery (slope on x ~ 0.8)
  expect_equal(unname(coef(fit)["x", 1]), 0.8, tolerance = 0.15)
  ## per-component sum-to-zero holds exactly
  expect_lt(abs(sum(fit$car$phi[, 1])), 1e-6)
})

test_that("G > 1 CAR mixture recovers regimes, intercepts, and phi", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 700, G = 2, tau = 0.5, spatial_error = TRUE,
                   lattice = 7, sep = 4, car_rho = 1.2, seed = 11)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 2, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 4L, seed = 2))
  cl <- apply(fit$posterior, 1L, which.max)
  acc <- max(mean(cl == d$truth$z), mean((3 - cl) == d$truth$z))
  expect_gt(acc, 0.75)
  expect_gt(diff(range(fit$beta_const[1, ])), 2)   # intercepts separated (true gap 4)
  ## phi recovered (allow a label swap)
  c_aligned <- max(
    mean(c(cor(fit$car$phi[, 1], d$truth$phi[, 1]),
           cor(fit$car$phi[, 2], d$truth$phi[, 2]))),
    mean(c(cor(fit$car$phi[, 1], d$truth$phi[, 2]),
           cor(fit$car$phi[, 2], d$truth$phi[, 1]))))
  expect_gt(c_aligned, 0.7)
})

test_that("phi_surface and predict broadcast the CAR effect", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 300, G = 1, tau = 0.5, spatial_error = TRUE,
                   lattice = 6, seed = 5)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  ps <- phi_surface(fit)
  expect_true(all(c("unit", "regime", "phi") %in% names(ps)))
  expect_equal(nrow(ps), nrow(fit$car$phi) * fit$G)
  qh <- predict(fit, type = "quantile")
  expect_length(qh, nrow(d$data))
})

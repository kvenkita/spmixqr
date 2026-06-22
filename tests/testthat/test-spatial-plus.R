test_that("Spatial+ reduces spatial-confounding bias under check loss", {
  d <- sim_spatial_confound(n = 400, beta = 1.0, confound = 0.9, seed = 11)
  Wn <- spq_weights(d$coords, type = "nngp", m = 10)
  ctl <- spmixqr_control(nstart = 1L, seed = 1)
  b <- function(f) coef(f)["x", 1]
  f_naive <- spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5,
                     spatial_error = FALSE, spatial_coef = FALSE, variance = "none", control = ctl)
  f_car <- spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5, spatial_error = TRUE,
                   spatial_coef = FALSE, spatial_W = Wn, variance = "none", control = ctl)
  f_sp <- suppressWarnings(spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5,
                  spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wn,
                  spatial_plus = TRUE, variance = "none", control = ctl))
  bias <- function(f) abs(b(f) - 1.0)
  expect_lt(bias(f_sp), bias(f_car))       # Spatial+ beats spatial-only ...
  expect_lt(bias(f_car), bias(f_naive))    # ... which beats naive
  expect_lt(bias(f_sp), 0.25)              # and lands near the truth
  r2 <- f_sp$diagnostics$spatial_plus$spatialR2
  expect_true(all(r2 >= 0 & r2 <= 1, na.rm = TRUE))
})

test_that("Spatial+ negative control: no confounding => slope undistorted", {
  d <- sim_spatial_confound(n = 350, beta = 1.0, confound = 0, seed = 5)
  Wn <- spq_weights(d$coords, type = "nngp", m = 10)
  ctl <- spmixqr_control(nstart = 1L, seed = 1)
  f_car <- spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5, spatial_error = TRUE,
                   spatial_coef = FALSE, spatial_W = Wn, variance = "none", control = ctl)
  f_sp <- spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5, spatial_error = TRUE,
                  spatial_coef = FALSE, spatial_W = Wn, spatial_plus = TRUE,
                  variance = "none", control = ctl)
  expect_lt(abs(coef(f_sp)["x", 1] - coef(f_car)["x", 1]), 0.1)
  expect_lt(abs(coef(f_sp)["x", 1] - 1.0), 0.15)
})

test_that("Spatial+ predict applies the stored smooths (predict == fitted)", {
  d <- sim_spatial_confound(n = 250, seed = 7)
  Wn <- spq_weights(d$coords, type = "nngp", m = 8)
  f_sp <- suppressWarnings(spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5,
                  spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wn,
                  spatial_plus = TRUE, variance = "none",
                  control = spmixqr_control(nstart = 1L, seed = 1)))
  pf <- predict(f_sp, newdata = d$data, newcoords = d$coords, type = "quantile")
  expect_lt(max(abs(as.numeric(pf) - f_sp$fitted_q)), 1e-8)
  ## predict must demand newcoords (else raw covariates vs deconfounded slopes)
  expect_error(predict(f_sp, newdata = d$data, type = "quantile"), "newcoords")
})

test_that("Spatial+ deconfounds via the resolution gap (FWL): rich smooth beats coarse", {
  ## Frisch-Waugh-Lovell: a residualisation smooth that matches the (penalised) spatial
  ## term is a no-op. Deconfounding therefore grows with the smooth's resolution. A coarse
  ## smooth (low k, ~ the spatial-error resolution) must remove LESS bias than a rich one.
  d <- sim_spatial_confound(n = 400, beta = 1.0, confound = 0.9, freq = 6, seed = 11)
  Wn <- spq_weights(d$coords, type = "nngp", m = 10)
  ctl <- spmixqr_control(nstart = 1L, seed = 1)
  bias <- function(k) {
    f <- suppressWarnings(spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wn,
                 spatial_plus = TRUE, spatial_plus_k = k, variance = "none", control = ctl))
    abs(coef(f)["x", 1] - 1.0)
  }
  expect_gt(bias(5L), bias(80L))     # coarse (no-op-ish) leaves more bias than rich
})

test_that("spatial_plus = FALSE leaves the fit unchanged (back-compat)", {
  d <- sim_spatial_confound(n = 200, seed = 3)
  ctl <- spmixqr_control(nstart = 1L, seed = 1)
  f0 <- spmixqr(y ~ x, d$data, coords = d$coords, G = 1, tau = 0.5,
                spatial_error = FALSE, spatial_coef = FALSE, variance = "none", control = ctl)
  expect_false(isTRUE(f0$spatial_plus))
  expect_null(f0$spatial_plus_smooths)
  expect_null(f0$diagnostics$spatial_plus)
})

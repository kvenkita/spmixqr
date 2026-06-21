## Point-data CAR coverage (the path that hid C1/M1/M2): build point coordinates by
## mapping the lattice CAR DGP's units onto a 2-D grid, fit with a knn weights graph
## via point coords, and exercise the default sandwich variance + out-of-sample predict.

## map a lattice CAR sim to per-observation point coordinates (r, c on the grid)
make_point_car <- function(n = 500L, G = 1L, lattice = 7L, seed = 7L, car_rho = 1.5) {
  d <- sim_spmixqr(n = n, G = G, tau = 0.5, spatial_error = TRUE,
                   lattice = lattice, car_rho = car_rho, seed = seed)
  side <- as.integer(lattice)
  u <- d$truth$unit                              # observation -> unit (1..L)
  ## expand.grid(r = 1:side, c = 1:side): r varies fastest, matching the sim
  rr <- ((u - 1L) %% side) + 1L
  cc <- ((u - 1L) %/% side) + 1L
  coords <- cbind(sx = as.numeric(rr), sy = as.numeric(cc))
  ## per-unit coordinates (1..L) in unit order, for aligning truth$phi
  ur <- ((seq_len(side * side) - 1L) %% side) + 1L
  uc <- ((seq_len(side * side) - 1L) %/% side) + 1L
  list(d = d, coords = coords, unit_coords = cbind(ur, uc), side = side)
}

test_that("point CAR (matrix coords + knn) fits, recovers phi, and predicts", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  pc <- make_point_car(n = 500L, G = 1L, lattice = 7L, seed = 7L)
  Wp <- spq_weights(pc$coords[!duplicated(pc$coords), , drop = FALSE],
                    type = "knn", k = 4L)
  fit <- spmixqr(y ~ x, pc$d$data, coords = pc$coords, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wp,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_true(isTRUE(fit$spatial_error))
  expect_identical(fit$car$units$mode, "point")

  ## phi recovery: align fitted units (by their stored coords) to the truth units
  fc <- fit$car$units$coords
  ord <- match(paste(pc$unit_coords[, 1], pc$unit_coords[, 2]),
               paste(fc[, 1], fc[, 2]))
  expect_false(anyNA(ord))
  phi_hat <- fit$car$phi[ord, 1]
  expect_gt(cor(phi_hat, pc$d$truth$phi[, 1]), 0.6)
  ## per-component sum-to-zero
  expect_lt(abs(sum(fit$car$phi[, 1])), 1e-6)

  ## predict at training coordinates (in-sample) works
  qh <- predict(fit, type = "quantile")
  expect_length(qh, nrow(pc$d$data))
  ## predict with newdata at exact training coords works (M1 fix)
  nd <- pc$d$data[1:5, , drop = FALSE]
  nc <- pc$coords[1:5, , drop = FALSE]
  qn <- predict(fit, newdata = nd, newcoords = nc, type = "quantile")
  expect_length(qn, 5L)
  expect_true(all(is.finite(qn)))
  ## predict at a slightly perturbed coordinate falls back to nearest training unit
  nc2 <- nc + 0.01
  qn2 <- predict(fit, newdata = nd, newcoords = nc2, type = "quantile")
  expect_true(all(is.finite(qn2)))
})

test_that("point CAR accepts coords = c('sx','sy') column names (M2 fix)", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  pc <- make_point_car(n = 400L, G = 1L, lattice = 6L, seed = 5L)
  dat <- cbind(pc$d$data, sx = pc$coords[, "sx"], sy = pc$coords[, "sy"])
  Wp <- spq_weights(pc$coords[!duplicated(pc$coords), , drop = FALSE],
                    type = "knn", k = 4L)
  fit <- spmixqr(y ~ x, dat, coords = c("sx", "sy"), G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wp,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_identical(fit$car$units$mode, "point")
  expect_equal(nrow(fit$car$phi), nrow(unique(pc$coords)))
})

test_that("point CAR with the DEFAULT variance='sandwich' returns finite SEs (C1 fix)", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  pc <- make_point_car(n = 400L, G = 1L, lattice = 6L, seed = 5L)
  Wp <- spq_weights(pc$coords[!duplicated(pc$coords), , drop = FALSE],
                    type = "knn", k = 4L)
  ## variance defaults to "sandwich" -- must not crash on the sparse CAR design
  fit <- spmixqr(y ~ x, pc$d$data, coords = pc$coords, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wp,
                 control = spmixqr_control(nstart = 1L, seed = 1))
  expect_identical(fit$se_method, "sandwich")
  V <- vcov(fit)
  expect_equal(V$method, "sandwich")
  ## beta SEs (the non-CAR coefficient block) are finite
  s <- summary(fit)
  expect_true(all(is.finite(s$comp[[1]][, "Std. Error"])))
  ci <- confint(fit)
  expect_true(all(is.finite(ci[[1]])))
})

test_that("NNGP precision is symmetric, PSD, and reduces to the dense Matern precision", {
  set.seed(1)
  co <- matrix(runif(40), 20, 2)
  nn <- spmixqr:::nngp_precision(co, m = 6, range = 0.3, nu = 0.5)
  expect_true(Matrix::isSymmetric(nn$Q))
  expect_true(Matrix::isSymmetric(nn$Wadj))
  expect_gt(min(eigen(as.matrix(nn$Q), only.values = TRUE)$values), -1e-8)  # PSD/PD
  ## dense limit: m = n-1 recovers solve(Matern covariance) (nu = 0.5, jitter-safe)
  co2 <- matrix(runif(24), 12, 2)
  Qf <- spmixqr:::nngp_precision(co2, m = 11, range = 0.4, nu = 0.5)$Q
  Cd <- spmixqr:::matern_cor(as.matrix(dist(co2)), 0.4, 0.5); diag(Cd) <- 1
  expect_lt(max(abs(as.matrix(Qf) - solve(Cd))), 1e-4)
})

test_that("NNGP sparsity is O(n*m^2) and grows ~linearly in n", {
  set.seed(2)
  n <- 60L; m <- 10L
  nn <- spmixqr:::nngp_precision(matrix(runif(2 * n), n, 2), m = m, range = 0.3)
  expect_lte(length(nn$Q@x), 4 * n * m)            # NOT (m+1)*n; the B'D^-1 B fill-in
  expect_gt(length(nn$Q@x), (m + 1) * n)           # genuinely denser than the edge bound
})

test_that("spq_weights(type='nngp') builds a usable proper-precision object", {
  set.seed(3)
  cc <- matrix(runif(80), 40, 2)
  w <- spq_weights(cc, type = "nngp", m = 6, nu = 0.5)
  expect_s3_class(w, "spq_weights")
  expect_identical(w$type, "nngp")
  expect_true(!is.null(w$Q) && Matrix::isSymmetric(w$Q))
  expect_equal(nrow(w$Q), 40)
  expect_output(print(w), "nngp")
})

test_that("NNGP folds into the penalised spatial-error M-step (unconstrained phi)", {
  data(meuse_zinc)
  cc <- as.matrix(meuse_zinc[, c("x", "y")])
  Wn <- spq_weights(cc, type = "nngp", m = 8)
  fit <- spmixqr(log(zinc) ~ dist, meuse_zinc, coords = cc, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = Wn,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_identical(fit$car$kind, "nngp")
  ph <- phi_surface(fit)
  expect_true(all(is.finite(ph$phi)))
  expect_equal(nrow(ph), nrow(cc))
  ## proper GP => phi is NOT forced to sum to zero (unlike ICAR)
  expect_gt(abs(sum(ph$phi)), 1e-6)
  expect_false(is.null(fit$diagnostics$moran))      # synthetic W drives Moran
})

test_that("non-NNGP weights and spatial_plus=FALSE are unchanged (back-compat)", {
  data(columbus_crime); data(columbus_W)
  W <- spq_weights(columbus_W, "supplied")
  reg <- factor(columbus_crime$id)
  f1 <- spmixqr(crime ~ income, columbus_crime, coords = reg, G = 1, tau = 0.5,
                spatial_error = TRUE, spatial_coef = FALSE, spatial_W = W,
                variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  expect_identical(f1$car$kind, "car")
  expect_false(isTRUE(f1$spatial_plus))
})

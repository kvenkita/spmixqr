test_that("Moran power: autocorrelation significant before, not after the CAR term", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 600, G = 1, tau = 0.5, spatial_error = TRUE,
                   lattice = 8, car_rho = 2, seed = 101)
  fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
  mo <- moran_resid(fit, nsim = 499)
  expect_s3_class(mo, "spq_moran")
  expect_lt(mo$before$p_value, 0.05)               # significant autocorr before
  expect_gt(mo$after$p_value, 0.05)                # absorbed after the CAR term
  expect_gt(mo$before$statistic, mo$after$statistic)
})

test_that("Moran size: no autocorrelation is near nominal (small Monte Carlo)", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  reps <- 20L; pvals <- numeric(reps)
  side <- 8L; L <- side * side
  cg <- expand.grid(a = seq_len(side), b = seq_len(side))
  W <- matrix(0, L, L)
  for (i in seq_len(L)) for (j in seq_len(L))
    if (i < j && abs(cg$a[i] - cg$a[j]) + abs(cg$b[i] - cg$b[j]) == 1) {
      W[i, j] <- 1; W[j, i] <- 1
    }
  spw <- spq_weights(W, type = "supplied", ids = as.character(seq_len(L)))
  for (r in seq_len(reps)) {
    set.seed(2000 + r)
    n <- 400; unit <- sample.int(L, n, replace = TRUE)
    x <- rnorm(n); y <- 1 + 0.5 * x + rnorm(n, 0, 0.8)  # NO spatial structure
    dat <- data.frame(y = y, x = x, region = as.character(unit))
    f <- suppressMessages(spmixqr(y ~ x, dat, coords = dat$region, G = 1, tau = 0.5,
           spatial_error = TRUE, spatial_coef = FALSE, spatial_W = spw,
           variance = "none", control = spmixqr_control(nstart = 1L, seed = 1)))
    pvals[r] <- moran_resid(f, nsim = 199)$before$p_value
  }
  ## rejection rate should be roughly nominal (loose band for a 20-rep MC)
  expect_lt(mean(pvals < 0.05), 0.25)
  expect_gt(mean(pvals), 0.2)
})

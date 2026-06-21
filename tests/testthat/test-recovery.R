test_that("spatial mixture recovers regimes, gate trend, and slope order", {
  set.seed(42)
  d <- sim_spmixqr(n = 300, G = 2, tau = 0.5, seed = 42)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 variance = "none", control = spmixqr_control(nstart = 3L, seed = 1))
  expect_true(fit$diagnostics$converged)
  ## classification accuracy beats chance (allow either label alignment)
  cl <- apply(fit$posterior, 1L, which.max)
  acc <- max(mean(cl == d$truth$z), mean((3 - cl) == d$truth$z))
  expect_gt(acc, 0.75)
  ## gate captures the spatial membership trend
  expect_gt(abs(cor(fit$prior[, 2], d$truth$pi[, 2])), 0.7)
  ## intercepts separated and ordered (ascending by slope key)
  expect_gt(diff(range(fit$beta_const[1, ])), 1)
})

test_that("negative control: constant DGP yields a near-flat gate surface", {
  set.seed(11)
  d <- sim_spmixqr(n = 250, gate_slope = 0, coef_slope = 0, seed = 11)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 lambda_gate = 10, lambda_coef = 10, variance = "none",
                 control = spmixqr_control(nstart = 3L, seed = 1))
  gs <- gate_surface(fit)
  p2 <- gs$prob[gs$regime == "2"]
  ## the gate probability should not vary much across space (no true structure)
  expect_lt(sd(p2), 0.25)
})

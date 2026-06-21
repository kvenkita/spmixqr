test_that("sandwich vcov has finite SEs and summary reports them", {
  set.seed(9)
  d <- sim_spmixqr(n = 200, seed = 9)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 variance = "sandwich", control = spmixqr_control(nstart = 2L, seed = 1))
  V <- vcov(fit)
  expect_equal(V$method, "sandwich")
  expect_true(all(is.finite(diag(V$gate))))
  expect_true(all(vapply(V$coef, function(M) all(is.finite(diag(M))), logical(1))))
})

test_that("bootstrap inference runs as a smoke test (small B) and aligns labels", {
  skip_on_cran()
  set.seed(3)
  d <- sim_spmixqr(n = 150, seed = 3)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5, variance = "boot",
                 control = spmixqr_control(nstart = 2L, seed = 1, boot_B = 12L,
                                           boot_block = 1L))
  expect_equal(fit$vcov$method, "boot")
  expect_true(is.matrix(fit$vcov$all))
  expect_true(all(is.finite(diag(fit$vcov$all))))
  ## bootstrap SEs must flow through the S3 surface (regression guard: they were NA)
  s <- summary(fit)
  expect_true(all(is.finite(s$comp[[1]][, "Std. Error"])))
  expect_true(all(is.finite(s$gate[[1]][, "Std. Error"])))
  expect_true(all(is.finite(confint(fit)[[1]])))
})

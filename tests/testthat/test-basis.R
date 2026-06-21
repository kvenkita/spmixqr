test_that("spatial basis builds, penalty is PSD, predicts at new locations", {
  set.seed(1)
  cc <- cbind(runif(120), runif(120))
  b <- spmixqr_basis(cc, type = "tp", k = 20)
  expect_s3_class(b, "spmixqr_basis")
  expect_equal(nrow(b$B), 120)
  expect_true(isSymmetric(unname(b$Omega)))
  ## PSD (not PD): smallest eigenvalue ~ 0, allowing tiny negative rounding
  ev <- eigen(b$Omega, only.values = TRUE)$values
  expect_gt(min(ev), -1e-6)
  expect_true(any(abs(ev) < 1e-6))
  Bn <- spmixqr:::predict_basis(b, rbind(c(0.5, 0.5), c(0.2, 0.8)))
  expect_equal(ncol(Bn), ncol(b$B))
})

test_that("k is clamped to unique locations", {
  set.seed(2)
  cc <- cbind(runif(8), runif(8))
  expect_message(b <- spmixqr_basis(cc, type = "tp", k = 20), "clamped")
  expect_lte(b$r, 7)
})

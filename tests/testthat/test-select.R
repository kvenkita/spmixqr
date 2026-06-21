test_that("spmixqr_select runs over a grid and returns a best fit (BIC)", {
  set.seed(1)
  d <- sim_spmixqr(n = 160, seed = 1)
  sel <- spmixqr_select(y ~ x, d$data, coords = d$coords, tau = 0.5, G_grid = 2,
                        lambda_gate_grid = c(1, 10), lambda_coef_grid = c(1, 10),
                        criterion = "bic", control = spmixqr_control(nstart = 2L, seed = 1))
  expect_s3_class(sel$fit, "spmixqr")
  expect_true(nrow(sel$table) == 4L)
  expect_true(is.finite(sel$best$score))
})

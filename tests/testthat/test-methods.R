test_that("S3 surface runs and predict types are coherent", {
  set.seed(9)
  d <- sim_spmixqr(n = 200, seed = 9)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 variance = "sandwich", control = spmixqr_control(nstart = 2L, seed = 1))
  expect_output(print(fit))
  s <- summary(fit); expect_s3_class(s, "summary.spmixqr")
  expect_output(print(s))
  expect_equal(dim(coef(fit)), c(2L, 2L))
  expect_true(is.list(vcov(fit)))
  ci <- confint(fit); expect_length(ci, 2L)
  expect_true(is.finite(AIC(fit)) && is.finite(BIC(fit)))
  expect_length(fitted(fit), 200L)
  expect_length(residuals(fit), 200L)
  expect_equal(nobs(fit), 200L)
  ## predict types
  pr <- predict(fit, type = "prob"); expect_equal(dim(pr), c(200L, 2L))
  cl <- predict(fit, type = "class"); expect_length(cl, 200L)
  qm <- predict(fit, type = "quantile"); expect_length(qm, 200L)
  qb <- predict(fit, type = "quantile_byclass"); expect_equal(dim(qb), c(200L, 2L))
  ## predict at new locations
  nd <- data.frame(x = rnorm(5)); ncoo <- cbind(runif(5), runif(5))
  expect_equal(nrow(predict(fit, nd, ncoo, type = "prob")), 5L)
})

test_that("coef_surface and gate_surface return tidy frames", {
  set.seed(10)
  d <- sim_spmixqr(n = 180, seed = 10)
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 variance = "none", control = spmixqr_control(nstart = 2L, seed = 1))
  cs <- coef_surface(fit); expect_true(all(c("regime", "slope") %in% names(cs)))
  gs <- gate_surface(fit); expect_true(all(c("regime", "prob") %in% names(gs)))
})

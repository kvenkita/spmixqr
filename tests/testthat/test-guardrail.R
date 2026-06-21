test_that("explicit spatial_gate = TRUE with spatial_error = TRUE raises a guardrail error", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 200, G = 2, tau = 0.5, spatial_error = TRUE, lattice = 5, seed = 3)
  expect_error(
    spmixqr(y ~ x, d$data, coords = d$region, G = 2, tau = 0.5, spatial_error = TRUE,
            spatial_gate = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
            variance = "none", control = spmixqr_control(nstart = 1L)),
    "not separately identified")
})

test_that("default spatial_gate flips to FALSE with a message under spatial_error", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 250, G = 2, tau = 0.5, spatial_error = TRUE, lattice = 5, seed = 4)
  expect_message(
    f <- spmixqr(y ~ x, d$data, coords = d$region, G = 2, tau = 0.5,
                 spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
                 variance = "none", control = spmixqr_control(nstart = 1L, seed = 1)),
    "turning the spatial gate off")
  expect_false(f$spatial_gate)
})

test_that("spmixqr_select threads spatial_error without tripping the guardrail", {
  skip_if_not_installed("spdep"); skip_if_not_installed("Matrix")
  d <- sim_spmixqr(n = 300, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 6, seed = 8)
  sel <- suppressMessages(spmixqr_select(
    y ~ x, d$data, coords = d$region, tau = 0.5, G_grid = 1,
    spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
    lambda_error_grid = c(0.5, 5), criterion = "bic",
    control = spmixqr_control(nstart = 1L, seed = 1)))
  expect_true(is.finite(sel$best$score))
  expect_true("lambda_error" %in% names(sel$best))
  expect_true(isTRUE(sel$fit$spatial_error))
})

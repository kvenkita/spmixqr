test_that("spq_weights builds symmetric W from supplied / queen / rook / distance / knn", {
  skip_if_not_installed("Matrix")
  ## supplied chain 1-2-3-4
  Wm <- matrix(0, 4, 4)
  Wm[1, 2] <- Wm[2, 1] <- Wm[2, 3] <- Wm[3, 2] <- Wm[3, 4] <- Wm[4, 3] <- 1
  w <- spq_weights(Wm, type = "supplied")
  expect_s3_class(w, "spq_weights")
  expect_true(Matrix::isSymmetric(w$W))
  expect_equal(w$n_comp, 1L)

  skip_if_not_installed("spdep")
  set.seed(1); pts <- cbind(runif(40), runif(40))
  wk <- spq_weights(pts, type = "knn", k = 4)
  expect_true(Matrix::isSymmetric(wk$W))           # knn directed -> symmetrised
  wd <- spq_weights(pts, type = "distance", d1 = 0, d2 = 0.4)
  expect_true(Matrix::isSymmetric(wd$W))

  skip_if_not_installed("sf")
  suppressMessages(library(sf))
  g <- st_make_grid(st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 3, ymax = 3))),
                    n = c(3, 3))
  gsf <- st_sf(id = 1:9, geometry = g)
  wq <- spq_weights(gsf, type = "queen")
  wr <- spq_weights(gsf, type = "rook")
  expect_gt(length(wq$W@x), length(wr$W@x))        # queen >= rook links
  expect_true(Matrix::isSymmetric(wq$W))
})

test_that("supplied W validators reject non-square / negative / asymmetric", {
  expect_error(spq_weights(matrix(1, 2, 3), type = "supplied"), "square")
  Wneg <- matrix(c(0, -1, -1, 0), 2, 2)
  expect_error(spq_weights(Wneg, type = "supplied"), "nonnegative")
  Wasym <- matrix(c(0, 1, 0, 0), 2, 2)
  expect_warning(spq_weights(Wasym, type = "supplied"), "asymmetric")
})

test_that("proper Q is PD; ICAR null dimension equals the number of components", {
  skip_if_not_installed("Matrix")
  Wm <- matrix(0, 4, 4)
  Wm[1, 2] <- Wm[2, 1] <- Wm[2, 3] <- Wm[3, 2] <- Wm[3, 4] <- Wm[4, 3] <- 1
  w <- spq_weights(Wm, type = "supplied")
  Q <- spmixqr:::make_car_precision(w, alpha = 0.95, car = "proper")
  ev <- eigen(as.matrix(Q), only.values = TRUE)$values
  expect_gt(min(ev), 0)                            # proper CAR is PD
  expect_equal(min(ev), 0.05, tolerance = 1e-6)    # min eig = 1 - alpha

  ## disconnected 2-block: ICAR null dim = 2 = number of components
  Wd <- matrix(0, 4, 4); Wd[1, 2] <- Wd[2, 1] <- 1; Wd[3, 4] <- Wd[4, 3] <- 1
  wd <- spq_weights(Wd, type = "supplied")
  expect_equal(wd$n_comp, 2L)
  Qi <- spmixqr:::make_car_precision(wd, car = "icar", eps = 0)
  null_dim <- sum(abs(eigen(as.matrix(Qi), only.values = TRUE)$values) < 1e-8)
  expect_equal(null_dim, 2L)
})

test_that("constraint absorption gives per-component sum-to-zero contrasts", {
  skip_if_not_installed("Matrix")
  Wm <- matrix(0, 4, 4)
  Wm[1, 2] <- Wm[2, 1] <- Wm[2, 3] <- Wm[3, 2] <- Wm[3, 4] <- Wm[4, 3] <- 1
  w <- spq_weights(Wm, type = "supplied")
  Q <- spmixqr:::make_car_precision(w, 0.95, "proper")
  R <- spmixqr:::incidence_matrix(c(1, 1, 2, 3, 4, 4), 4)
  mem <- spmixqr:::components_from_W(w$W)$membership
  ab <- spmixqr:::absorb_car_constraint(R, Q, mem)
  expect_equal(ab$Lp, 3L)                          # L - n_comp = 4 - 1
  expect_true(all(abs(colSums(as.matrix(ab$Tmat))) < 1e-8))
})

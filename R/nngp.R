## Nearest-Neighbour Gaussian Process (NNGP) sparse Matern precision for the
## point-referenced spatial-error term. Builds the Vecchia/NNGP precision
## Q = (I - B)' D^{-1} (I - B) (Datta, Banerjee, Finley & Gelfand 2016, JASA), a sparse,
## proper, full-rank GMRF precision that plugs into the same penalised spatial-error
## M-step used by the CAR path. Pure spdep + Matrix (no spNNGP/INLA dependency).
##
## Design decisions verified in the adversarial spec review (2026-06-22):
##  - PROPER full-rank precision => phi is fit UNCONSTRAINED (no ICAR sum-to-zero).
##  - max-min ordering (Guinness 2018) is the default; ordering affects the
##    approximation but not PSD, and m -> n-1 recovers the dense Matern precision.
##  - default range = a fraction of the domain extent (NOT a median-NN multiple, which
##    shrinks with n); range is in RAW coordinate units. range/nu are fixed-or-selected
##    (weakly identified from one realisation; Zhang 2004), like lambda_phi.
##  - relative jitter (~1e-10 * mean(diag)) + a conditional-variance floor keep the
##    neighbour solves well-conditioned (Matern-1.5 is ill-conditioned).

#' Matern correlation (closed forms for nu = 0.5 or 1.5).
#' @keywords internal
matern_cor <- function(h, range, nu = 0.5) {
  h <- h / range
  if (abs(nu - 0.5) < 1e-8) exp(-h)
  else if (abs(nu - 1.5) < 1e-8) (1 + sqrt(3) * h) * exp(-sqrt(3) * h)
  else stop("NNGP supports nu in {0.5, 1.5}.", call. = FALSE)
}

#' Max-min ordering of points (Guinness 2018): greedily pick the point farthest
#' (in min-distance) from those already chosen. O(n^2) pure R; fine to a few thousand.
#' @keywords internal
maxmin_order <- function(coords) {
  n <- nrow(coords)
  if (n <= 2L) return(seq_len(n))
  ctr <- colMeans(coords)
  d0 <- sqrt(rowSums((coords - rep(ctr, each = n))^2))
  ord <- integer(n)
  ord[1] <- which.min(d0)                       # start nearest the centroid
  mind <- sqrt(colSums((t(coords) - coords[ord[1], ])^2))
  mind[ord[1]] <- -Inf
  for (i in 2:n) {
    pick <- which.max(mind)
    ord[i] <- pick
    di <- sqrt(colSums((t(coords) - coords[pick, ])^2))
    mind <- pmin(mind, di)
    mind[pick] <- -Inf
  }
  ord
}

#' NNGP sparse Matern precision and its neighbour graph.
#'
#' Builds `Q = (I - B)' D^{-1} (I - B)` (Datta et al. 2016) on the supplied point
#' coordinates, returning the sparse precision (in the input row order) plus the
#' symmetrised neighbour adjacency (used by diagnostics / the spatial-block bootstrap).
#'
#' @param coords An `n x 2` coordinate matrix (distinct locations; raw units).
#' @param m Number of nearest earlier-neighbours (default 10).
#' @param range Matern range in raw coordinate units; `NULL` => ~0.1 * domain extent.
#' @param nu Matern smoothness in `{0.5, 1.5}`.
#' @return A list with `Q` (sparse symmetric `dgCMatrix`, input order), `Wadj` (sparse
#'   symmetric 0/1 neighbour adjacency), `range`, `m`, `nu`, `ordering` (max-min order).
#' @references Datta, Banerjee, Finley & Gelfand (2016, JASA); Guinness (2018,
#'   Technometrics); Zhang (2004, JASA).
#' @keywords internal
nngp_precision <- function(coords, m = 10L, range = NULL, nu = 0.5) {
  coords <- as.matrix(coords); storage.mode(coords) <- "double"
  n <- nrow(coords); m <- as.integer(min(m, n - 1L))
  if (n < 2L) stop("NNGP needs at least 2 distinct locations.", call. = FALSE)
  if (is.null(range)) {
    ext <- max(stats::dist(coords[seq_len(min(n, 500L)), , drop = FALSE]))  # domain extent
    range <- 0.1 * ext
  }
  if (!is.finite(range) || range <= 0) stop("NNGP `range` must be positive.", call. = FALSE)
  ord <- maxmin_order(coords)
  co  <- coords[ord, , drop = FALSE]
  Bi <- integer(0); Bj <- integer(0); Bx <- numeric(0)
  Ai <- integer(0); Aj <- integer(0)                 # neighbour graph (ordered indices)
  D  <- numeric(n); c11 <- 1
  D[1] <- c11
  for (i in 2:n) {
    prev <- seq_len(i - 1L)
    di   <- sqrt(colSums((t(co[prev, , drop = FALSE]) - co[i, ])^2))
    nb   <- prev[order(di)[seq_len(min(m, i - 1L))]]
    hnn  <- as.matrix(stats::dist(co[nb, , drop = FALSE]))
    Cnn  <- matrix(matern_cor(hnn, range, nu), length(nb), length(nb)); diag(Cnn) <- c11
    cin  <- matern_cor(sqrt(colSums((t(co[nb, , drop = FALSE]) - co[i, ])^2)), range, nu)
    jit  <- 1e-10 * mean(diag(Cnn))
    bi   <- solve(Cnn + jit * diag(length(nb)), cin)
    D[i] <- max(c11 - sum(bi * cin), 1e-8 * c11)     # conditional variance floor
    Bi <- c(Bi, rep(i, length(nb))); Bj <- c(Bj, nb); Bx <- c(Bx, bi)
    Ai <- c(Ai, rep(i, length(nb))); Aj <- c(Aj, nb)
  }
  ImB  <- Matrix::Diagonal(n) - Matrix::sparseMatrix(i = Bi, j = Bj, x = Bx, dims = c(n, n))
  Qord <- Matrix::crossprod(ImB, Matrix::Diagonal(x = 1 / D) %*% ImB)
  inv  <- order(ord)                                  # ordered -> input index
  Q    <- Matrix::forceSymmetric(Qord[inv, inv])
  Q    <- methods::as(methods::as(Q, "CsparseMatrix"), "generalMatrix")
  ## symmetric 0/1 neighbour adjacency in INPUT order (for Moran / bootstrap blocks)
  if (length(Ai)) {
    oi <- ord[Ai]; oj <- ord[Aj]                      # map ordered -> input indices
    Wadj <- Matrix::sparseMatrix(i = c(oi, oj), j = c(oj, oi),
                                 x = rep(1, 2 * length(oi)), dims = c(n, n))
    Wadj@x[] <- 1                                     # collapse any double-counts to 1
  } else Wadj <- Matrix::Matrix(0, n, n, sparse = TRUE)
  Wadj <- methods::as(methods::as(Matrix::forceSymmetric(Wadj), "CsparseMatrix"),
                      "generalMatrix")
  list(Q = Q, Wadj = Wadj, range = range, m = m, nu = nu, ordering = ord)
}

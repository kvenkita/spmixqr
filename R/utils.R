## Small internal helpers. Some are vendored (re-implemented) from the mixqr /
## mixqrgate internals, which are not exported; credited at each site.

#' Numeric floor for densities / probabilities (avoid log(0)).
#' @keywords internal
.dens_floor <- 1e-12

#' Row-normalise a non-negative matrix so rows sum to one.
#' @keywords internal
normalize_rows <- function(M) {
  rs <- rowSums(M)
  rs[rs <= 0] <- 1
  M / rs
}

#' Residual matrix (response minus fitted quantiles), one column per regime.
#' `fitted` is an n by G matrix of fitted quantiles.
#' @keywords internal
residual_matrix_fitted <- function(y, fitted) y - fitted

#' Asymmetric-Laplace density at the tau-th quantile parametrisation
#' (Yu & Moyeed 2001), scale `sigma`. Vendored from the mixqr ALD idiom.
#' @keywords internal
ald_density <- function(u, sigma, tau) {
  sigma <- max(sigma, 1e-8)
  (tau * (1 - tau) / sigma) * exp(-rho_tau(u, tau) / sigma)
}

#' Check (pinball) loss: u times (tau minus the indicator that u is negative).
#' @keywords internal
rho_tau <- function(u, tau) u * (tau - (u < 0))

#' Standard-normal cdf/pdf shorthands.
#' @keywords internal
.pnorm <- function(x) stats::pnorm(x)
#' @keywords internal
.dnorm <- function(x) stats::dnorm(x)

#' Symmetrise a matrix (numerical hygiene for covariances).
#'
#' Coerces a \pkg{Matrix} S4 object to a base dense matrix first, so the result is
#' always a plain numeric matrix (covariances downstream are read with base `diag()`).
#' @keywords internal
symmetrise <- function(A) {
  if (isS4(A) || methods::is(A, "Matrix")) A <- as.matrix(A)
  (A + t(A)) / 2
}

#' Minimal Moore-Penrose pseudo-inverse (avoids a MASS dependency).
#' Vendored from mixqrgate's MASS_ginv.
#' @keywords internal
ginv_small <- function(A, tol = 1e-10) {
  s <- svd(A)
  d <- s$d
  dinv <- ifelse(d > tol * max(d), 1 / d, 0)
  s$v %*% (dinv * t(s$u))
}

#' Safe solve with a ridge fallback for near-singular systems.
#' @keywords internal
safe_solve <- function(A, b, ridge = 1e-8) {
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (is.null(out)) out <- solve(A + ridge * diag(nrow(A)), b)
  out
}

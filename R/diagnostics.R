## Spatial diagnostics for the CAR-error module. Permutation Moran's I on the
## responsibility-weighted quantile residual aggregated to the spatial unit, reported
## before and after the CAR term. Permutation inference (not Gaussian-ML asymptotics)
## is robust under the ALD working likelihood.

#' Permutation Moran's I on the spatial-unit residuals of an spmixqr fit
#'
#' Computes Moran's I of the responsibility-weighted quantile residual aggregated to
#' the spatial unit, with a permutation p-value (random relabelling of units). For a
#' CAR fit (`spatial_error = TRUE`) it reports the residual Moran's I *after* the CAR
#' term and, by re-fitting without it, the value *before* - the negative-control /
#' power diagnostic. For a non-spatial-error fit it reports the single residual Moran's
#' I against the supplied weights.
#'
#' The mixture residual is `r_s = sum_{i in s} sum_k p_ik (y_i - x_i'beta_k - phi_{k,s})
#' / n_s`, the responsibility-weighted residual aggregated to unit `s` (per-regime
#' residuals are ambiguous off-support; documented in the spec).
#'
#' @param object A fitted [spmixqr] object.
#' @param spatial_W A [spq_weights()] object or weights matrix. Defaults to the fit's
#'   CAR weights (`object$car$W`) when present; otherwise required.
#' @param nsim Number of permutations for the p-value.
#' @return A list of class `spq_moran` with `before`/`after` (each a list with
#'   `statistic`, `p_value`, `n_units`), or a single `statistic`/`p_value` for a
#'   non-CAR fit.
#' @references \pkg{spdep} permutation Moran's I (Bivand et al.); Cliff & Ord (1981).
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE,
#'                  lattice = 6, seed = 1)
#' fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
#'                spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
#'                variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
#' moran_resid(fit, nsim = 199)
#' }
#' @export
moran_resid <- function(object, spatial_W = NULL, nsim = 999L) {
  if (!inherits(object, "spmixqr")) stop("`object` must be an spmixqr fit.")
  spqw <- resolve_moran_weights(object, spatial_W)
  unit_idx <- moran_unit_index(object, spqw)

  ## residual AFTER (the current fit, including any CAR term)
  r_after <- unit_residual(object, unit_idx, nrow(spqw$W))
  after <- perm_moran(r_after, spqw$W, nsim)

  if (!isTRUE(object$spatial_error)) {
    out <- list(statistic = after$statistic, p_value = after$p_value,
                n_units = after$n_units, nsim = nsim)
    class(out) <- "spq_moran"
    return(out)
  }

  ## residual BEFORE: the same fit's residual with the CAR effect removed (phi = 0),
  ## i.e. the residual against the beta surface alone. This isolates the spatial
  ## autocorrelation the CAR term absorbs without a fragile call re-evaluation.
  before <- tryCatch({
    car_block <- object$design$car_block
    beta_only <- if (length(car_block) > 0L)
      object$coefficients[-car_block, , drop = FALSE] else object$coefficients
    Xt <- object$design$Xt
    if (length(car_block) > 0L) Xt <- Xt[, -car_block, drop = FALSE]
    fq0 <- as.matrix(Xt %*% beta_only)               # n x G, no phi
    r_before <- unit_residual_internal(object$y, fq0, object$posterior,
                                       unit_idx, nrow(spqw$W))
    perm_moran(r_before, spqw$W, nsim)
  }, error = function(e) NULL)

  out <- list(before = before, after = after, nsim = nsim,
              note = "permutation Moran's I on responsibility-weighted unit residuals")
  class(out) <- "spq_moran"
  out
}

#' @keywords internal
resolve_moran_weights <- function(object, spatial_W) {
  if (!is.null(spatial_W)) {
    if (inherits(spatial_W, "spq_weights")) return(spatial_W)
    return(spq_weights(spatial_W, type = "supplied"))
  }
  if (!is.null(object$car)) return(object$car$W)
  stop("Supply `spatial_W` (the fit has no CAR weights to reuse).", call. = FALSE)
}

#' Observation -> unit index for the Moran aggregation.
#' @keywords internal
moran_unit_index <- function(object, spqw) {
  if (!is.null(object$car)) return(object$car$units$unit_idx)
  ## non-CAR fit: derive units from coords/areal matching spqw dimension
  geo <- object$coords
  L <- nrow(spqw$W)
  if (!is.null(geo) && identical(geo$mode, "areal")) {
    region <- as.factor(geo$region)
    return(as.integer(region))
  }
  if (!is.null(geo) && identical(geo$mode, "point")) {
    cc <- geo$coords
    key <- apply(cc, 1L, function(z) paste(z, collapse = "_"))
    return(match(key, unique(key)))
  }
  ## fallback: each row its own unit (requires L == n)
  if (L == length(object$y)) return(seq_along(object$y))
  stop("Cannot map observations to the supplied weights' units.", call. = FALSE)
}

#' Responsibility-weighted residual aggregated to the spatial unit (length L).
#' @keywords internal
unit_residual <- function(object, unit_idx, L) {
  unit_residual_internal(object$y, object$fitted_q, object$posterior, unit_idx, L)
}

#' Array-level unit residual (used during the fit before the object exists).
#' @keywords internal
unit_residual_internal <- function(y, fitted_q, posterior, unit_idx, L) {
  e <- y - fitted_q                        # n x G residuals (fitted_q includes phi)
  rw <- rowSums(posterior * e)             # responsibility-weighted residual, length n
  f <- factor(unit_idx, levels = seq_len(L))
  num <- as.numeric(tapply(rw, f, sum))
  cnt <- as.numeric(tapply(rep(1, length(rw)), f, sum))
  num[is.na(num)] <- 0; cnt[is.na(cnt)] <- 0
  ifelse(cnt > 0, num / pmax(cnt, 1), NA_real_)
}

#' Permutation Moran's I for a unit-level vector against a weights matrix.
#' @keywords internal
perm_moran <- function(x, W, nsim = 999L) {
  ok <- is.finite(x)
  if (sum(ok) < 3L)
    return(list(statistic = NA_real_, p_value = NA_real_, n_units = sum(ok)))
  W <- W[ok, ok, drop = FALSE]
  x <- x[ok]
  n <- length(x)
  S0 <- sum(W)
  if (S0 <= 0)
    return(list(statistic = NA_real_, p_value = NA_real_, n_units = n))
  moran_I <- function(z) {
    zc <- z - mean(z)
    num <- as.numeric(Matrix::crossprod(zc, W %*% zc))
    den <- sum(zc^2)
    (n / S0) * (num / den)
  }
  obs <- moran_I(x)
  perm <- vapply(seq_len(nsim), function(i) moran_I(x[sample.int(n)]), numeric(1))
  ## two-sided permutation p-value (rank of |obs| among permutations)
  p <- (1 + sum(abs(perm) >= abs(obs))) / (nsim + 1)
  list(statistic = obs, p_value = p, n_units = n)
}

#' @export
print.spq_moran <- function(x, ...) {
  cat("Permutation Moran's I (responsibility-weighted unit residuals)\n")
  fmt <- function(s) sprintf("I = %.4f   p = %.4f   (units = %d)",
                             s$statistic, s$p_value, s$n_units)
  if (!is.null(x$before) || !is.null(x$after)) {
    if (!is.null(x$before)) cat("  before CAR: ", fmt(x$before), "\n", sep = "")
    if (!is.null(x$after))  cat("  after  CAR: ", fmt(x$after), "\n", sep = "")
  } else {
    cat("  ", fmt(x), "\n", sep = "")
  }
  invisible(x)
}

## Spatial+ confounding safeguard (Dupont, Wood & Augustin 2022, Biometrics).
## Residualise each covariate against a spatial smooth and fit the quantile model on the
## residuals, so a smoothly-spatial covariate no longer competes with the spatial random
## effect (phi). It is loss-agnostic (operates only on covariates) so it composes with the
## quantile EM and with the CAR / NNGP spatial-error term.
##
## IMPORTANT (Frisch-Waugh-Lovell, verified in spec review): residualising on a basis that
## matches the model's spatial term is a NO-OP. Spatial+ deconfounds only to the extent the
## residualisation smooth out-resolves the (penalised/shrunk) spatial term. The default
## residualisation basis is therefore made richer than the typical spatial-error resolution.
## The residualisation is mean (least-squares) based; under check loss it deconfounds the
## tau-slope cleanly for symmetric error (validated empirically) -- the reported slope is the
## effect of the NON-SPATIAL part of the covariate.

#' Residualise covariates against a spatial smooth (Spatial+ stage 1).
#'
#' For each non-intercept column of the component design, regress it on a spatial smooth
#' (a thin-plate smooth of point coordinates, or a Markov-random-field smooth over areal
#' regions) and replace it by the residual. Returns the residualised design, the fitted
#' smooths (needed to residualise `newdata` at prediction), and the per-covariate spatial
#' R-squared removed.
#'
#' @param X The component design matrix (with intercept).
#' @param slope_idx Integer indices of the non-intercept columns to residualise.
#' @param geo The resolved geography (`resolve_coords()` output): point `coords` or areal
#'   `region` + `areal` nb.
#' @param k Basis dimension for the smooth (`NULL` = a generous default richer than the
#'   spatial-error resolution).
#' @param w Optional numeric vector of length `nrow(X)` of observation weights;
#'   `NULL` (default) reproduces the unweighted residualization.
#' @return A list with `X` (residualised design), `smooths` (named list of fitted `gam`s,
#'   one per residualised covariate; `NULL` for skipped constant columns), and `spatialR2`
#'   (a data frame of the spatial R-squared removed per covariate).
#' @references Dupont, Wood & Augustin (2022, Biometrics).
#' @keywords internal
spatial_plus_residualize <- function(X, slope_idx, geo, k = NULL, w = NULL) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Package 'mgcv' is required for spatial_plus.", call. = FALSE)
  nm <- colnames(X)[slope_idx]
  smooths <- stats::setNames(vector("list", length(slope_idx)), nm)
  r2 <- rep(NA_real_, length(slope_idx))
  point <- geo$mode == "point"
  if (point) {
    cc <- as.matrix(geo$coords)
    base_df <- data.frame(c1 = cc[, 1], c2 = cc[, 2])
    ## richer-than-spatial-error default (FWL: must out-resolve the penalised phi)
    kk <- if (is.null(k)) min(100L, max(30L, nrow(cc) %/% 4L)) else as.integer(k)
  } else {
    if (is.null(geo$areal))
      stop("Spatial+ for areal data needs `areal` (an spdep nb).", call. = FALSE)
    reg <- factor(as.character(geo$region))
    nbL <- lapply(geo$areal, as.integer)
    rid <- attr(geo$areal, "region.id")
    names(nbL) <- if (!is.null(rid)) as.character(rid) else levels(reg)
    base_df <- data.frame(reg = reg)
  }
  for (jj in seq_along(slope_idx)) {
    j <- slope_idx[jj]; xj <- as.numeric(X[, j])
    if (stats::var(xj) < .Machine$double.eps) next        # constant covariate: skip
    df <- base_df; df$xj <- xj
    if (!is.null(w)) df$.w <- as.numeric(w)
    g <- tryCatch(
      if (point) mgcv::gam(xj ~ s(c1, c2, bs = "tp", k = kk), data = df, method = "REML",
                           weights = if (!is.null(w)) df$.w else NULL)
      else       mgcv::gam(xj ~ s(reg, bs = "mrf", xt = list(nb = nbL)), data = df,
                           method = "REML", weights = if (!is.null(w)) df$.w else NULL),
      error = function(e) NULL)
    if (is.null(g)) next
    res <- xj - as.numeric(stats::fitted(g))
    r2[jj] <- 1 - stats::var(res) / stats::var(xj)
    X[, j] <- res
    smooths[[jj]] <- g
  }
  list(X = X, smooths = smooths,
       spatialR2 = data.frame(variable = nm, spatialR2 = round(r2, 4),
                              row.names = NULL))
}

#' Apply stored Spatial+ smooths to residualise new covariates (predict stage).
#'
#' Mirrors [spatial_plus_residualize()] on `newdata`: subtracts each stored smooth's
#' prediction (at the new coordinates / regions) from the corresponding covariate column.
#' @keywords internal
spatial_plus_apply <- function(Xnew, smooths, geo_new) {
  if (length(smooths) == 0L) return(Xnew)
  point <- geo_new$mode == "point"
  if (point) {
    cc <- as.matrix(geo_new$coords); nd <- data.frame(c1 = cc[, 1], c2 = cc[, 2])
  } else {
    nd <- data.frame(reg = factor(as.character(geo_new$region)))
  }
  for (nm in names(smooths)) {
    g <- smooths[[nm]]
    if (is.null(g) || !(nm %in% colnames(Xnew))) next
    Xnew[, nm] <- as.numeric(Xnew[, nm]) - as.numeric(stats::predict(g, newdata = nd))
  }
  Xnew
}

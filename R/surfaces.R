#' Spatial coefficient surfaces
#'
#' Evaluates each regime's covariate-effect slope surface(s) `beta_kj(s)` over a set
#' of locations (the training locations by default), returning a tidy data frame for
#' mapping.
#'
#' @param object A fitted [spmixqr] object.
#' @param newcoords Optional coordinates (point) or region labels (areal); defaults
#'   to the training locations.
#' @param covariate Which slope covariate (index among the non-intercept terms).
#' @return A data frame with coordinates, `regime`, and `slope`.
#' @export
coef_surface <- function(object, newcoords = NULL, covariate = 1L) {
  if (!object$spatial_coef)
    stop("This model has flat (non-spatial) component slopes.")
  if (is.null(newcoords)) {
    B <- object$design$Xt[, object$design$spat_blocks[[covariate]], drop = FALSE]
    ## recover the raw basis B from the augmented columns is x_j*B; instead use stored basis
    B <- object$basis$B
    cc <- if (object$coords$mode == "point") object$coords$coords else object$coords$region
  } else {
    loc <- if (object$basis$type == "mrf") newcoords else as.matrix(newcoords)
    B <- predict_basis(object$basis, loc)
    cc <- newcoords
  }
  const_pos <- object$design$const_rows[1L + covariate]
  blk <- object$design$spat_blocks[[covariate]]
  G <- object$G
  vals <- lapply(seq_len(G), function(k) {
    s <- object$coefficients[const_pos, k] + as.numeric(B %*% object$coefficients[blk, k])
    data.frame(regime = factor(k), slope = s)
  })
  out <- do.call(rbind, Map(function(df, k) df, vals, seq_len(G)))
  if (object$basis$type != "mrf") {
    ccm <- as.matrix(cc)
    out$coord1 <- rep(ccm[, 1L], G); out$coord2 <- rep(ccm[, 2L], G)
  } else {
    out$region <- rep(as.character(cc), G)
  }
  out
}

#' CAR spatial-error surfaces
#'
#' Returns the per-regime CAR spatial random effect `phi_k` (the mean-zero
#' spatial-level surface) as a tidy data frame for mapping, mirroring [coef_surface()]
#' and [gate_surface()]. One row per (unit, regime).
#'
#' @param object A fitted [spmixqr] object with `spatial_error = TRUE`.
#' @param newunits Optional unit identifiers to restrict / reorder the output
#'   (validated against the CAR unit ids). `NULL` uses all training units.
#' @return A data frame with `unit`, `regime`, and `phi`.
#' @examples
#' \donttest{
#' d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 6,
#'                  seed = 1)
#' fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
#'                spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
#'                variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
#' head(phi_surface(fit))
#' }
#' @export
phi_surface <- function(object, newunits = NULL) {
  if (!isTRUE(object$spatial_error) || is.null(object$car))
    stop("This model has no CAR spatial-error term (`spatial_error = FALSE`).")
  phi <- object$car$phi
  ids <- object$car$units$ids
  if (!is.null(newunits)) {
    sel <- match(as.character(newunits), ids)
    if (anyNA(sel)) stop("`newunits` contains units not seen in training.")
    phi <- phi[sel, , drop = FALSE]; ids <- ids[sel]
  }
  G <- ncol(phi); L <- nrow(phi)
  data.frame(unit = factor(rep(ids, times = G), levels = ids),
             regime = factor(rep(seq_len(G), each = L)),
             phi = as.numeric(phi))
}

#' Spatial gate surfaces
#'
#' Evaluates the mixing probabilities over space, returning a tidy data frame.
#'
#' @param object A fitted [spmixqr] object.
#' @param newcoords Optional coordinates / region labels; defaults to training.
#' @param newdata Optional data frame for gating covariates (if any).
#' @return A data frame with coordinates, `regime`, and `prob`.
#' @export
gate_surface <- function(object, newcoords = NULL, newdata = NULL) {
  if (is.null(newcoords)) {
    Z <- object$design$Z
    cc <- if (object$coords$mode == "point") object$coords$coords else object$coords$region
  } else {
    if (is.null(newdata) && !identical(object$gating, ~1))
      stop("This model has gating covariates; supply `newdata` with them for new-location gate surfaces.")
    W <- if (is.null(newdata)) matrix(1, NROW(newcoords), 1)
         else stats::model.matrix(object$gating, stats::model.frame(object$gating, newdata))
    loc <- if (object$basis$type == "mrf") newcoords else as.matrix(newcoords)
    B <- if (object$spatial_gate) predict_basis(object$basis, loc) else NULL
    Z <- if (object$spatial_gate) cbind(W, B) else W
    cc <- newcoords
  }
  prob <- gate_predict(object$gamma, Z)
  G <- object$G; n <- nrow(prob)
  out <- data.frame(regime = factor(rep(seq_len(G), each = n)),
                    prob = as.numeric(prob))
  if (object$basis$type != "mrf") {
    ccm <- as.matrix(cc)
    out$coord1 <- rep(ccm[, 1L], G); out$coord2 <- rep(ccm[, 2L], G)
  } else {
    out$region <- rep(as.character(cc), G)
  }
  out
}

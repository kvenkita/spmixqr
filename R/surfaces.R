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

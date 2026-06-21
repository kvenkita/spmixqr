#' Predictions from a spatial quantile mixture model
#'
#' @param object A fitted [spmixqr] object.
#' @param newdata Optional data frame of predictors. If `NULL`, in-sample.
#' @param newcoords Coordinates (point: two-column) or region labels (areal) for
#'   `newdata`. Required when `newdata` is supplied and the model is spatial.
#' @param type One of `"prob"` (gate probabilities over space, `n x G`), `"class"`
#'   (most probable regime by the gate), `"posterior"` (responsibilities; needs the
#'   response in `newdata`), `"quantile"` (gate-mixed conditional quantile), or
#'   `"quantile_byclass"` (per-regime conditional quantile, `n x G`).
#' @param ... Ignored.
#' @return A vector or matrix depending on `type`.
#' @export
predict.spmixqr <- function(object, newdata = NULL, newcoords = NULL,
                            type = c("prob", "class", "posterior", "quantile",
                                     "quantile_byclass"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    Xt <- object$design$Xt; Z <- object$design$Z
    prob <- object$prior; post <- object$posterior
    fq <- object$fitted_q
  } else {
    d <- predict_designs(object, newdata, newcoords)
    Xt <- d$Xt; Z <- d$Z
    prob <- gate_predict(object$gamma, Z)
    fq <- Xt %*% object$coefficients
    post <- NULL
    if (type == "posterior") {
      yv <- tryCatch(stats::model.response(stats::model.frame(object$formula, newdata)),
                     error = function(e) NULL)
      if (is.null(yv)) stop("type='posterior' needs the response in `newdata`.")
      dm <- density_matrix(object, yv - fq)
      post <- normalize_rows(pmax(prob * dm, .dens_floor))
    }
  }
  switch(type,
    prob = prob,
    class = apply(prob, 1L, which.max),
    posterior = post,
    quantile = rowSums(prob * fq),
    quantile_byclass = fq)
}

#' Build component / gate designs at new data (same column layout as the fit).
#' @keywords internal
predict_designs <- function(object, newdata, newcoords) {
  tt <- stats::delete.response(object$terms)
  X <- stats::model.matrix(tt, stats::model.frame(tt, newdata))
  W <- stats::model.matrix(object$gating, stats::model.frame(object$gating, newdata))
  B <- NULL
  if (object$spatial_gate || object$spatial_coef) {
    if (is.null(newcoords)) stop("Spatial model: supply `newcoords` for prediction.")
    loc <- if (object$basis$type == "mrf") newcoords else as.matrix(newcoords)
    B <- predict_basis(object$basis, loc)
  }
  des <- build_designs(X, W, B, object$basis$Omega, object$design$slope_idx,
                       object$spatial_gate, object$spatial_coef,
                       object$lambda_gate, object$lambda_coef,
                       object$gate_ridge %||% 1e-3,
                       object$design$r)
  list(Xt = des$Xt, Z = des$Z)
}

#' Per-regime component densities at residual matrix `e` (n x G).
#' @keywords internal
density_matrix <- function(object, e) {
  G <- object$G; n <- nrow(e)
  if (object$method == "ald" && !is.null(object$sigma)) {
    vapply(seq_len(G), function(k) ald_density(e[, k], object$sigma[k], object$tau),
           numeric(n))
  } else if (!is.null(object$dens)) {
    vapply(seq_len(G), function(k) pmax(object$dens[[k]]$eval(e[, k]), .dens_floor),
           numeric(n))
  } else {
    vapply(seq_len(G), function(k) ald_density(e[, k], 1, object$tau), numeric(n))
  }
}

#' Null-coalescing helper.
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

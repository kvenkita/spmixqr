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
    ## beta-only fitted quantile (CAR columns excluded; add the phi surface below)
    beta_only <- if (length(object$design$car_block) > 0L)
      object$coefficients[-object$design$car_block, , drop = FALSE]
    else object$coefficients
    fq <- as.matrix(Xt %*% beta_only)
    if (isTRUE(object$spatial_error)) fq <- fq + predict_car_phi(object, newcoords)
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
  ## Spatial+: residualise new covariates with the stored smooths (else predictions
  ## silently mix raw covariates with deconfounded slopes -> wrong quantiles).
  if (isTRUE(object$spatial_plus) && length(object$spatial_plus_smooths)) {
    if (is.null(newcoords))
      stop("Spatial+ model: supply `newcoords` so covariates can be residualised.",
           call. = FALSE)
    geo_new <- list(mode = object$coords$mode,
                    coords = if (object$coords$mode == "point") newcoords else NULL,
                    region = if (object$coords$mode == "areal") newcoords else NULL)
    X <- spatial_plus_apply(X, object$spatial_plus_smooths, geo_new)
  }
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

#' Broadcast the CAR phi surface (L x G) to new observations (n x G) via their unit.
#' @keywords internal
predict_car_phi <- function(object, newcoords) {
  car <- object$car
  if (is.null(newcoords))
    stop("Spatial-error model: supply `newcoords` (region labels or coordinates).",
         call. = FALSE)
  ids <- car$units$ids
  if (identical(car$units$mode, "areal")) {
    u <- match(as.character(newcoords), ids)
    if (anyNA(u))
      stop("predict for a spatial-error model needs `newcoords` matching the training ",
           "spatial units (CAR phi is defined on the fixed unit set).", call. = FALSE)
  } else {
    cc <- as.matrix(newcoords)
    storage.mode(cc) <- "double"
    tc <- car$units$coords
    if (is.null(tc))
      stop("predict for a point spatial-error model needs stored training unit ",
           "coordinates (refit with spmixqr 0.2.0 or later).", call. = FALSE)
    ## CAR phi is defined on the fixed training unit set. Match each new point to a
    ## training unit: exact coordinate match first, else the nearest training-unit
    ## centroid by Euclidean distance (documented nearest-unit broadcast).
    tkey <- apply(tc, 1L, function(z) paste(z, collapse = "_"))
    nkey <- apply(cc, 1L, function(z) paste(z, collapse = "_"))
    u <- match(nkey, tkey)
    miss <- which(is.na(u))
    if (length(miss)) {
      for (i in miss) {
        dvec <- sqrt(rowSums((tc - matrix(cc[i, ], nrow(tc), ncol(tc),
                                          byrow = TRUE))^2))
        u[i] <- which.min(dvec)
      }
    }
  }
  car$phi[u, , drop = FALSE]
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

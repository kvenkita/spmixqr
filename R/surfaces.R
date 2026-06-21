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

#' CAR spatial-error surfaces (with uncertainty)
#'
#' Returns the per-regime CAR spatial random effect `phi_k` (the mean-zero
#' spatial-level deviation, on the response scale at the fitted quantile) as a tidy
#' data frame for mapping, mirroring [coef_surface()] and [gate_surface()]. One row per
#' (unit, regime). With `ci = TRUE` it adds standard errors and a confidence interval,
#' and flags the `credible` units whose interval excludes zero (reliable hot/cold
#' spots, as distinct from smoothing). With `scale = "exp"` it adds the multiplicative
#' deviation `mult = exp(phi)` (interpretable for log-outcome models: e.g. `mult` 1.18
#' means the outcome runs about 18% above what the covariates predict there).
#'
#' @param object A fitted [spmixqr] object with `spatial_error = TRUE`.
#' @param newunits Optional unit identifiers to restrict / reorder the output
#'   (validated against the CAR unit ids). `NULL` uses all training units.
#' @param ci Add `se`, `lower`, `upper`, and `credible` columns. Requires the fit to
#'   carry a covariance (`variance = "boot"`, recommended, or `"sandwich"`).
#' @param level Confidence level for the interval.
#' @param scale `"link"` (the response/quantile scale, default) or `"exp"` (add the
#'   multiplicative deviation `exp(phi)`, for log-outcome models).
#' @return A data frame with `unit`, `regime`, `phi`, and (per `ci`/`scale`) `se`,
#'   `lower`, `upper`, `credible`, `mult`, `mult_lower`, `mult_upper`.
#' @examples
#' \donttest{
#' d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 6,
#'                  seed = 1)
#' fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
#'                spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
#'                control = spmixqr_control(nstart = 1L, seed = 1))
#' head(phi_surface(fit, ci = TRUE))
#' }
#' @export
phi_surface <- function(object, newunits = NULL, ci = FALSE, level = 0.95,
                        scale = c("link", "exp")) {
  if (!isTRUE(object$spatial_error) || is.null(object$car))
    stop("This model has no CAR spatial-error term (`spatial_error = FALSE`).")
  scale <- match.arg(scale)
  phi <- object$car$phi
  ids <- object$car$units$ids
  se <- if (ci) phi_se(object) else NULL
  sel <- if (is.null(newunits)) seq_along(ids) else match(as.character(newunits), ids)
  if (anyNA(sel)) stop("`newunits` contains units not seen in training.")
  phi <- phi[sel, , drop = FALSE]; ids <- ids[sel]
  if (!is.null(se)) se <- se[sel, , drop = FALSE]
  G <- ncol(phi); L <- nrow(phi)
  out <- data.frame(unit = factor(rep(ids, times = G), levels = ids),
                    regime = factor(rep(seq_len(G), each = L)),
                    phi = as.numeric(phi))
  if (ci) {
    if (is.null(se)) {
      warning("No standard errors available (fit with variance = \"boot\" or ",
              "\"sandwich\"); CI columns are NA.", call. = FALSE)
      out$se <- NA_real_; out$lower <- NA_real_; out$upper <- NA_real_
    } else {
      z <- stats::qnorm(1 - (1 - level) / 2)
      out$se <- as.numeric(se)
      out$lower <- out$phi - z * out$se
      out$upper <- out$phi + z * out$se
      ## "credible": the interval excludes zero (a reliable hot/cold spot)
      out$credible <- is.finite(out$se) & (out$lower > 0 | out$upper < 0)
    }
  }
  if (scale == "exp") {
    out$mult <- exp(out$phi)             # multiplicative deviation (log-outcome models)
    if (!is.null(out$lower)) { out$mult_lower <- exp(out$lower); out$mult_upper <- exp(out$upper) }
  }
  out
}

#' Standard errors of the CAR spatial effect phi (L x G).
#'
#' Propagates the fit's stored coefficient covariance through the sum-to-zero
#' constraint transform: `Var(phi_k) = T V_red,k T'`, with `V_red,k` the CAR sub-block
#' of the regime-k coefficient covariance and `T` the constraint-absorption basis. The
#' **bootstrap** (`variance = "boot"`) is recommended: it refits the penalised pipeline
#' and so reflects the shrinkage of the random effect, whereas the classification-
#' conditional sandwich is a fast alternative that ignores the penalty and can disagree.
#' With few spatial blocks the bootstrap intervals can be optimistic (small-sample
#' caveat). Returns `NULL` if no covariance is stored.
#' @keywords internal
phi_se <- function(object) {
  V <- object$vcov
  if (is.null(V) || is.null(V$coef)) return(NULL)
  Tm <- object$car$Tmat; cb <- object$car$car_block
  L <- nrow(object$car$phi); G <- ncol(object$car$phi)
  if (is.null(Tm) || length(cb) == 0L) return(matrix(0, L, G))
  Tm <- as.matrix(Tm)
  vapply(seq_len(G), function(k) {
    Vk <- V$coef[[k]]
    if (is.null(Vk)) return(rep(NA_real_, L))
    Vred <- as.matrix(Vk)[cb, cb, drop = FALSE]
    Vphi <- Tm %*% Vred %*% t(Tm)
    sqrt(pmax(diag(Vphi), 0))
  }, numeric(L))
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

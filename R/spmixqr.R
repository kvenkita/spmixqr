#' Fit a spatial finite mixture of quantile regressions
#'
#' Fits `G` latent regimes, each a `tau`-th quantile regression, whose mixing
#' probabilities vary over space through a spatial basis and (optionally) whose
#' covariate-effect slopes are spatially varying surfaces. Estimation is a penalised
#' EM. With `G = 1` this is a penalised spatially-varying-coefficient quantile
#' regression. See the package primer for a worked example.
#'
#' @param formula Component model formula, e.g. `y ~ x1 + x2`.
#' @param data A data frame.
#' @param coords Point coordinates: a two-column matrix/data frame, or the names of
#'   two columns in `data`. For areal data, the region label per row (a factor or a
#'   column name), paired with `areal`.
#' @param areal An \pkg{spdep} `nb`/`listw` neighbour object over the regions (areal
#'   data). `NULL` for point data.
#' @param G Number of regimes.
#' @param tau Quantile level in (0, 1).
#' @param gating Gating-covariate formula (besides space), e.g. `~ z`. Default `~ 1`.
#' @param spatial_gate Logical; let the mixing probabilities vary over space.
#' @param spatial_coef Logical; let component slopes vary over space (scalar
#'   intercepts; see Details).
#' @param method Component density: `"ald"` (asymmetric Laplace) or `"kde"` (Wu & Yao
#'   constrained kernel density).
#' @param lambda_gate,lambda_coef Roughness penalties for the gate and slope
#'   surfaces. `NULL` uses the control default; choose with [spmixqr_select()].
#' @param variance Inference: `"sandwich"` (fast, classification-conditional;
#'   default), `"boot"` (spatial-block/xy bootstrap; recommended for reporting), or
#'   `"none"`.
#' @param basis Advanced override: a prebuilt [spmixqr_basis()]. Usually `NULL`
#'   (built internally from `coords`/`areal`).
#' @param control A [spmixqr_control()] list.
#'
#' @details A free spatial intercept surface is not identified separately from the
#'   spatial gate (both move the marginal quantile level across space). v1 therefore
#'   models spatial level/membership through the gate and spatial covariate effects
#'   through the slope surfaces, with scalar per-regime intercepts. Component
#'   identities are ordered by a sign-stable separating slope; the `label_stability`
#'   diagnostic warns when slope surfaces cross (a single global label is then only
#'   approximately coherent).
#'
#' @return An object of class `spmixqr`.
#' @references Reich, Fuentes & Dunson (2011); Wu & Yao (2016); Fernandes, Guerre &
#'   Horta (2021).
#' @examples
#' set.seed(1)
#' d <- sim_spmixqr(n = 200, G = 2, tau = 0.5)
#' fit <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
#'                variance = "none", control = spmixqr_control(nstart = 2L))
#' fit
#' @export
spmixqr <- function(formula, data, coords = NULL, areal = NULL, G = 2L, tau = 0.5,
                    gating = ~1, spatial_gate = TRUE, spatial_coef = TRUE,
                    method = c("ald", "kde"),
                    lambda_gate = NULL, lambda_coef = NULL,
                    variance = c("sandwich", "boot", "none"),
                    basis = NULL, control = spmixqr_control()) {
  method <- match.arg(method)
  if (is.character(variance) || is.null(variance))
    variance <- if (is.null(variance)) "none" else match.arg(variance)
  if (!is.null(control$seed)) set.seed(control$seed)
  if (tau <= 0 || tau >= 1) stop("`tau` must be in (0, 1).")
  G <- as.integer(G)

  ## ---- response + component design ----
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  n <- length(y); p <- ncol(X)
  slope_idx <- which(colnames(X) != "(Intercept)")
  if (length(slope_idx) == 0L && spatial_coef)
    stop("`spatial_coef = TRUE` needs at least one non-intercept covariate.")

  ## ---- gating design ----
  Wmf <- stats::model.frame(gating, data)
  W <- stats::model.matrix(attr(Wmf, "terms"), Wmf)

  ## ---- spatial basis (built once, shared by gate and slope surfaces) ----
  need_basis <- spatial_gate || spatial_coef
  geo <- if (need_basis) resolve_coords(coords, areal, data, n) else NULL
  if (need_basis && is.null(basis)) {
    bt <- if (!is.null(areal)) "mrf" else control$basis_type
    locarg <- if (geo$mode == "areal") geo$region else geo$coords
    basis <- spmixqr_basis(locarg, areal = geo$areal, type = bt,
                           k = control$k, scale_coords = control$scale_coords)
  }
  B <- if (need_basis) basis$B else NULL
  Omega <- if (need_basis) basis$Omega else NULL
  r <- if (need_basis) basis$r else 0L
  lam_g <- if (is.null(lambda_gate)) 1 else lambda_gate
  lam_b <- if (is.null(lambda_coef)) 1 else lambda_coef

  ## ---- assemble augmented designs + penalties + column bookkeeping ----
  des <- build_designs(X, W, B, Omega, slope_idx, spatial_gate, spatial_coef,
                       lam_g, lam_b, control$gate_ridge, r)
  Xt <- des$Xt; Z <- des$Z; Pen_beta <- des$Pen_beta; Pen_gamma <- des$Pen_gamma

  kdectrl <- mixqr::mixqr_control(bandwidth = control$bandwidth,
                                  kde_grid = control$kde_grid)
  h_rate <- smooth_bandwidth(n, control$sm_scale)

  ## ---- multi-start EM ----
  coords_num <- if (need_basis && geo$mode == "point") geo$coords else NULL
  starts <- init_starts(y, X, coords_num, G, tau, control$nstart)
  fits <- lapply(starts, function(p0)
    tryCatch(spatial_em_fit(y, Xt, Z, G, tau, method, p0, Pen_beta, Pen_gamma,
                            h_rate, spatial_coef, kdectrl, control),
             error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0L) stop("All EM starts failed.")
  best <- fits[[which.max(vapply(fits, function(f) f$objective, numeric(1)))]]

  ## ---- relabel by sign-stable separating slope, then refit the gate ----
  beta_const <- extract_const(best$beta, des, p)
  ord <- order_components(beta_const, control$label_order, control$order_var)
  best <- permute_fit(best, ord)
  gfit <- pen_irls_multinom(Z, best$posterior, Pen_gamma, control$gate_maxit,
                            control$gate_tol)
  best$gamma <- gfit$gamma; best$prior <- gfit$pi; best$gate_hessian <- gfit$hessian
  beta_const <- extract_const(best$beta, des, p)

  ## ---- label-stability diagnostic (slope-surface crossing) ----
  lab_stab <- NA_real_
  if (spatial_coef && G > 1L && length(slope_idx) > 0L) {
    ss <- slope_surface_matrix(best$beta, des, B, control$order_var)
    lab_stab <- label_stability(ss)
    if (lab_stab > 0.2)
      warning(sprintf(paste("slope surfaces cross at %.0f%% of locations:",
                            "a single global label is only approximately coherent",
                            "(see ?spmixqr Details)."), 100 * lab_stab))
  }

  ## ---- effective df, AIC/BIC ----
  edf <- compute_edf(best, Xt, Z, Pen_beta, Pen_gamma, G, des, spatial_coef)
  ll <- best$loglik
  aic <- if (is.finite(ll)) -2 * ll + 2 * edf$total else NA_real_
  bic <- if (is.finite(ll)) -2 * ll + log(n) * edf$total else NA_real_

  ## ---- diagnostics ----
  occ <- colSums(best$posterior)
  ent <- mean(apply(best$posterior, 1L, function(pp)
    -sum(pp[pp > 0] * log(pp[pp > 0])) / log(G)))
  diagnostics <- list(converged = best$converged, n_starts = length(fits),
                      iters = best$iters, gate_cond = best$gate_cond,
                      occupancy = occ, class_entropy = if (G > 1L) ent else 0,
                      label_stability = lab_stab, spatial_edf = edf$spatial,
                      smoothing_h = best$h_used)

  obj <- structure(list(
    coefficients = best$beta, beta_const = beta_const, gamma = best$gamma,
    prior = best$prior, posterior = best$posterior, sigma = best$sigma,
    dens = best$dens, fitted_q = X_fitted(Xt, best$beta), residuals_m = best$dm,
    loglik = ll, edf = edf$total, edf_detail = edf, aic = aic, bic = bic,
    lambda_gate = lam_g, lambda_coef = lam_b, basis = basis, design = des,
    control = control, gate_ridge = control$gate_ridge,
    se_method = variance, vcov = NULL, diagnostics = diagnostics,
    call = match.call(), terms = attr(mf, "terms"), formula = formula,
    gating = gating, tau = tau, G = G, method = method,
    spatial_gate = spatial_gate, spatial_coef = spatial_coef,
    coords = if (need_basis) geo else NULL, X = X, W = W, y = y,
    h = best$h_used),
    class = "spmixqr")

  ## ---- inference ----
  if (variance == "sandwich") obj$vcov <- sandwich_vcov(obj)
  else if (variance == "boot") obj$vcov <- bootstrap_vcov(obj, data)
  obj
}

#' Resolve coordinates / areal specification.
#' @keywords internal
resolve_coords <- function(coords, areal, data, n) {
  if (!is.null(areal)) {
    region <- if (is.character(coords) && length(coords) == 1L) data[[coords]] else coords
    if (is.null(region)) stop("Areal fit needs `coords` = region labels (length n).")
    return(list(mode = "areal", region = as.factor(region), areal = areal,
                coords = NULL))
  }
  if (is.null(coords)) stop("Spatial terms need `coords` (point) or `areal` (areal).")
  cc <- if (is.character(coords)) as.matrix(data[, coords, drop = FALSE]) else as.matrix(coords)
  if (ncol(cc) != 2L) stop("`coords` must have two columns for point data.")
  storage.mode(cc) <- "double"
  list(mode = "point", coords = cc, areal = NULL, region = NULL)
}

#' Assemble augmented component / gate designs and their penalties.
#' @keywords internal
build_designs <- function(X, W, B, Omega, slope_idx, spatial_gate, spatial_coef,
                          lam_g, lam_b, gate_ridge, r) {
  n <- nrow(X)
  ## component design Xt and Pen_beta
  if (spatial_coef && length(slope_idx) > 0L && !is.null(B)) {
    cols <- list(X[, 1L, drop = FALSE]); const_rows <- 1L; spat_blocks <- list()
    pos <- 1L
    for (j in slope_idx) {
      cols[[length(cols) + 1L]] <- X[, j, drop = FALSE]   # constant slope
      pos <- pos + 1L; const_rows <- c(const_rows, pos)
      cols[[length(cols) + 1L]] <- X[, j] * B             # spatial deviation
      spat_blocks[[length(spat_blocks) + 1L]] <- (pos + 1L):(pos + r)
      pos <- pos + r
    }
    Xt <- do.call(cbind, cols); colnames(Xt) <- NULL
    P <- ncol(Xt)
    Pen_beta <- matrix(0, P, P)
    for (blk in spat_blocks) Pen_beta[blk, blk] <- lam_b * Omega
    intercept_row <- 1L
  } else {
    Xt <- X; P <- ncol(Xt)
    Pen_beta <- matrix(0, P, P)
    const_rows <- seq_len(P); spat_blocks <- list(); intercept_row <- 1L
  }
  ## gate design Z and Pen_gamma
  if (spatial_gate && !is.null(B)) {
    Z <- cbind(W, B); q1 <- ncol(Z)
    Pen_gamma <- matrix(0, q1, q1)
    qw <- ncol(W)
    Pen_gamma[seq_len(qw), seq_len(qw)] <- gate_ridge * diag(qw)
    gblk <- (qw + 1L):q1
    Pen_gamma[gblk, gblk] <- lam_g * Omega + gate_ridge * diag(r)  # ICAR null guard
    gate_spat <- gblk
  } else {
    Z <- W; q1 <- ncol(Z)
    Pen_gamma <- gate_ridge * diag(q1); gate_spat <- integer(0)
  }
  list(Xt = Xt, Z = Z, Pen_beta = Pen_beta, Pen_gamma = Pen_gamma,
       const_rows = const_rows, spat_blocks = spat_blocks,
       intercept_row = intercept_row, slope_idx = slope_idx,
       gate_spat = gate_spat, qw = ncol(W), r = r, spatial_coef = spatial_coef)
}

#' Constant coefficient matrix (p x G) from the augmented coefficients.
#' @keywords internal
extract_const <- function(beta, des, p) {
  cr <- des$const_rows
  if (length(cr) == nrow(beta)) return(beta)
  beta[cr, , drop = FALSE]
}

#' Permute a fit's components by an ordering.
#' @keywords internal
permute_fit <- function(fit, ord) {
  fit$beta <- fit$beta[, ord, drop = FALSE]
  fit$posterior <- fit$posterior[, ord, drop = FALSE]
  fit$prior <- fit$prior[, ord, drop = FALSE]
  fit$dm <- fit$dm[, ord, drop = FALSE]
  if (!is.null(fit$sigma)) fit$sigma <- fit$sigma[ord]
  if (!is.null(fit$dens)) fit$dens <- fit$dens[ord]
  fit
}

#' Ordering-covariate slope surface (n x G): beta_const + B' theta per regime.
#' @keywords internal
slope_surface_matrix <- function(beta, des, B, order_var) {
  j <- order_var
  if (j > length(des$spat_blocks)) j <- 1L
  const_pos <- des$const_rows[1L + j]
  blk <- des$spat_blocks[[j]]
  G <- ncol(beta)
  vapply(seq_len(G), function(k) beta[const_pos, k] + as.numeric(B %*% beta[blk, k]),
         numeric(nrow(B)))
}

#' Fitted conditional quantile per regime (n x G).
#' @keywords internal
X_fitted <- function(Xt, beta) Xt %*% beta

#' Effective degrees of freedom (component smoothers + gate smoother).
#' @keywords internal
compute_edf <- function(fit, Xt, Z, Pen_beta, Pen_gamma, G, des, spatial_coef) {
  comp <- 0
  for (k in seq_len(G)) {
    e <- fit$posterior[, k]                          # weight ~ responsibility
    A <- crossprod(Xt, e * Xt)
    H <- A + Pen_beta + 1e-8 * diag(ncol(Xt))
    comp <- comp + sum(diag(safe_solve(H, A)))
  }
  ## spatial portion = component edf beyond the constant (intercept + flat-slope) part
  n_const <- if (length(des$spat_blocks) > 0L) G * length(des$const_rows) else 0
  spatial <- max(0, comp - n_const)
  ## gate smoother edf
  gate_edf <- 0
  if (ncol(fit$gamma) > 0) {
    Ag <- -fit$gate_hessian
    ## unpenalised part = Ag - blkdiag(Pen_gamma)
    K <- ncol(fit$gamma); q1 <- ncol(Z)
    Pblk <- matrix(0, q1 * K, q1 * K)
    for (a in seq_len(K)) {
      ia <- ((a - 1L) * q1 + 1L):(a * q1); Pblk[ia, ia] <- Pen_gamma
    }
    gate_edf <- tryCatch(sum(diag(safe_solve(Ag, Ag - Pblk))), error = function(e) K)
  }
  list(total = comp + gate_edf + G, component = comp, gate = gate_edf,
       spatial = spatial)
}

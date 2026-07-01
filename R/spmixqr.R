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
#' @param spatial_error Logical; add a per-regime conditional-autoregressive (CAR)
#'   spatial random effect `phi_k` (a mean-zero spatial-level surface) to each
#'   component. When `TRUE` the spatial gate is turned off by default (a free spatial
#'   gate aliases with the CAR intercept surface; see Details); an explicit
#'   `spatial_gate = TRUE` then raises a guardrail error.
#' @param spatial_W A [spq_weights()] object or a square weights matrix for the CAR
#'   term. `NULL` auto-builds: queen contiguity from `areal` (an `nb`/`listw`), else a
#'   k-nearest-neighbour graph from point `coords`.
#' @param car CAR precision family: `"proper"` (Leroux `alpha(D-W)+(1-alpha)I`,
#'   default) or `"icar"` (intrinsic `D-W` with a sum-to-zero constraint).
#' @param car_alpha Proper-CAR spatial-dependence strength in `[0, 1]` (default
#'   `0.95`). Fixed by default (weakly identified; disclosed in `summary()`), not
#'   selected by the check loss.
#' @param spatial_plus Logical; apply the **Spatial+** confounding safeguard (Dupont,
#'   Wood & Augustin 2022): residualise each covariate against a spatial smooth and fit
#'   on the residuals, so a smoothly-spatial covariate no longer competes with the
#'   spatial random effect. The reported slopes are then the effect of the *non-spatial*
#'   part of each covariate. Deconfounds only to the extent the residualisation smooth
#'   out-resolves the (penalised) spatial term (Frisch--Waugh--Lovell); the default
#'   smooth is made richer for this reason. Composes with `spatial_error` (CAR or NNGP).
#' @param spatial_plus_k Basis dimension for the Spatial+ covariate smooth (`NULL` uses
#'   a generous default richer than the spatial-error resolution).
#' @param method Component density: `"ald"` (asymmetric Laplace) or `"kde"` (Wu & Yao
#'   constrained kernel density).
#' @param lambda_gate,lambda_coef,lambda_error Roughness / CAR penalties for the gate,
#'   slope surfaces, and CAR effect (`lambda_error` = `lambda_phi`). `NULL` uses the
#'   control default; choose with [spmixqr_select()].
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
                    spatial_error = FALSE, spatial_W = NULL,
                    car = c("proper", "icar"), car_alpha = 0.95,
                    spatial_plus = FALSE, spatial_plus_k = NULL,
                    method = c("ald", "kde"),
                    lambda_gate = NULL, lambda_coef = NULL, lambda_error = NULL,
                    variance = c("sandwich", "boot", "none"),
                    weights = NULL, weights_type = c("sampling", "frequency", "precision"),
                    basis = NULL, control = spmixqr_control()) {
  method <- match.arg(method)
  car <- match.arg(car)
  weights_type <- match.arg(weights_type)
  if (is.character(variance) || is.null(variance))
    variance <- if (is.null(variance)) "none" else match.arg(variance)
  if (!is.null(control$seed)) set.seed(control$seed)
  if (tau <= 0 || tau >= 1) stop("`tau` must be in (0, 1).")
  G <- as.integer(G)

  ## ---- spatial-error guardrail: a free spatial gate aliases with the CAR level ----
  if (isTRUE(spatial_error)) {
    if (!missing(spatial_gate) && isTRUE(spatial_gate))
      stop("`spatial_error = TRUE` cannot be combined with an explicit ",
           "`spatial_gate = TRUE`: a free spatial gate and the per-regime CAR ",
           "intercept surface are not separately identified (they both shift the ",
           "marginal quantile level across space). Use a covariate or constant gate ",
           "(`gating = ~1` or `gating = ~z`) with the CAR error, or drop the CAR ",
           "error. See ?spmixqr Details.", call. = FALSE)
    if (missing(spatial_gate) && isTRUE(spatial_gate)) {
      spatial_gate <- FALSE
      message("spatial_error = TRUE: turning the spatial gate off (the free spatial ",
              "gate aliases with the CAR intercept surface). Set gating = ~z for a ",
              "covariate gate.")
    }
  }

  ## ---- response + component design ----
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  n <- length(y); p <- ncol(X)
  wobj <- resolve_weights(weights, weights_type, data, n)
  w_obs <- wobj$w          # normalized (mean = 1); stored on fit and used for fitting
  w_fit <- w_obs
  slope_idx <- which(colnames(X) != "(Intercept)")
  if (length(slope_idx) == 0L && spatial_coef)
    stop("`spatial_coef = TRUE` needs at least one non-intercept covariate.")

  ## ---- gating design ----
  Wmf <- stats::model.frame(gating, data)
  W <- stats::model.matrix(attr(Wmf, "terms"), Wmf)

  ## ---- spatial basis (built once, shared by gate and slope surfaces) ----
  need_basis <- spatial_gate || spatial_coef
  need_geo <- need_basis || isTRUE(spatial_error) || isTRUE(spatial_plus)
  geo <- if (need_geo) resolve_coords(coords, areal, data, n,
                                      spatial_error = isTRUE(spatial_error) || isTRUE(spatial_plus),
                                      spatial_W = spatial_W) else NULL

  ## ---- Spatial+ confounding safeguard: residualise covariates against a spatial smooth ----
  sp_plus <- NULL
  if (isTRUE(spatial_plus)) {
    if (length(slope_idx) == 0L)
      stop("`spatial_plus = TRUE` needs at least one non-intercept covariate.", call. = FALSE)
    if (isTRUE(spatial_coef))
      warning("spatial_plus with spatial_coef = TRUE: the slope SURFACE then multiplies the ",
              "residualised (deconfounded) covariate, so its estimand also changes. v1 ",
              "applies Spatial+ to the constant-slope interpretation; read surfaces with care.",
              call. = FALSE)
    sp_plus <- spatial_plus_residualize(X, slope_idx, geo, k = spatial_plus_k)
    X <- sp_plus$X
  }
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
  lam_phi <- if (is.null(lambda_error)) control$lambda_error_default else lambda_error

  ## ---- CAR spatial-error structure (constraint-absorbed incidence + precision) ----
  car_obj <- NULL; car_list <- NULL
  if (isTRUE(spatial_error)) {
    car_obj <- resolve_car(spatial_W, geo, car, car_alpha, n)
    car_list <- list(Rt = car_obj$Rt, Qt = car_obj$Qt, Lp = car_obj$Lp,
                     lambda = lam_phi)
  }

  ## ---- assemble augmented designs + penalties + column bookkeeping ----
  des <- build_designs(X, W, B, Omega, slope_idx, spatial_gate, spatial_coef,
                       lam_g, lam_b, control$gate_ridge, r, car = car_list)
  Xt <- des$Xt; Z <- des$Z; Pen_beta <- des$Pen_beta; Pen_gamma <- des$Pen_gamma

  kdectrl <- mixqr::mixqr_control(bandwidth = control$bandwidth,
                                  kde_grid = control$kde_grid)
  h_rate <- smooth_bandwidth(n, control$sm_scale)

  ## ---- multi-start EM ----
  coords_num <- if (need_basis && geo$mode == "point") geo$coords else NULL
  starts <- init_starts(y, X, coords_num, G, tau, control$nstart)
  fits <- lapply(starts, function(p0)
    tryCatch(spatial_em_fit(y, Xt, Z, G, tau, method, p0, Pen_beta, Pen_gamma,
                            h_rate, spatial_coef, kdectrl, control,
                            spatial_error = isTRUE(spatial_error), w_obs = w_fit),
             error = function(e) NULL))
  fits <- Filter(Negate(is.null), fits)
  if (length(fits) == 0L) stop("All EM starts failed.")
  best <- fits[[which.max(vapply(fits, function(f) f$objective, numeric(1)))]]

  ## ---- relabel by sign-stable separating slope, then refit the gate ----
  beta_const <- extract_const(best$beta, des, p)
  ord <- order_components(beta_const, control$label_order, control$order_var)
  best <- permute_fit(best, ord)
  gfit <- pen_irls_multinom(Z, best$posterior, Pen_gamma, control$gate_maxit,
                            control$gate_tol, w = w_fit)
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
  edf <- compute_edf(best, Xt, Z, Pen_beta, Pen_gamma, G, des, spatial_coef, w_obs = w_fit)
  ll <- best$loglik
  aic <- if (is.finite(ll)) -2 * ll + 2 * edf$total else NA_real_
  bic <- if (is.finite(ll)) -2 * ll + log(n) * edf$total else NA_real_

  ## ---- CAR spatial-error effects (recover phi on the L units, mean-zero) ----
  car_slot <- NULL
  if (isTRUE(spatial_error)) {
    car_slot <- build_car_slot(best$beta, des$car_block, car_obj, car, car_alpha,
                               lam_phi)
  }

  ## ---- diagnostics ----
  occ <- colSums(best$posterior)
  ent <- mean(apply(best$posterior, 1L, function(pp)
    -sum(pp[pp > 0] * log(pp[pp > 0])) / log(G)))
  diagnostics <- list(converged = best$converged, n_starts = length(fits),
                      iters = best$iters, gate_cond = best$gate_cond,
                      occupancy = occ, class_entropy = if (G > 1L) ent else 0,
                      label_stability = lab_stab, spatial_edf = edf$spatial,
                      smoothing_h = best$h_used,
                      spatial_plus = if (!is.null(sp_plus)) sp_plus$spatialR2 else NULL)
  ## residual Moran's I (after the CAR term) for the summary block
  if (isTRUE(spatial_error) && !is.null(car_obj)) {
    diagnostics$moran <- tryCatch({
      r <- unit_residual_internal(y, X_fitted(Xt, best$beta), best$posterior,
                                  car_obj$unit_idx, nrow(car_obj$Q))
      perm_moran(r, car_obj$spqw$W, nsim = 199L)
    }, error = function(e) NULL)
  }

  obj <- structure(list(
    coefficients = best$beta, beta_const = beta_const, gamma = best$gamma,
    prior = best$prior, posterior = best$posterior, sigma = best$sigma,
    dens = best$dens, fitted_q = X_fitted(Xt, best$beta), residuals_m = best$dm,
    loglik = ll, edf = edf$total, edf_detail = edf, aic = aic, bic = bic,
    lambda_gate = lam_g, lambda_coef = lam_b, lambda_error = lam_phi,
    basis = basis, design = des,
    control = control, gate_ridge = control$gate_ridge,
    se_method = variance, vcov = NULL, diagnostics = diagnostics,
    call = match.call(), terms = attr(mf, "terms"), formula = formula,
    gating = gating, tau = tau, G = G, method = method,
    spatial_gate = spatial_gate, spatial_coef = spatial_coef,
    spatial_error = isTRUE(spatial_error), car = car_slot,
    spatial_plus = isTRUE(spatial_plus),
    spatial_plus_smooths = if (!is.null(sp_plus)) sp_plus$smooths else NULL,
    coords = if (need_geo) geo else NULL, X = X, W = W, y = y,
    weights = w_obs, prior_weights = wobj$raw, weights_type = weights_type,
    weights_sum = wobj$sum_raw,
    h = best$h_used),
    class = "spmixqr")

  ## ---- inference ----
  ## NNGP phi is unconstrained (L' = n units), so the classification-conditional sandwich
  ## densifies to O(n^2)/O(n^3); the spatial-block bootstrap is the scalable, recommended
  ## path for point-GP fits (synthesis item 7).
  if (variance == "sandwich" && !is.null(car_slot) && identical(car_slot$kind, "nngp") &&
      nrow(car_slot$phi) > 300L)
    warning("NNGP spatial-error fit with variance = 'sandwich' forms a dense ", nrow(car_slot$phi),
            "-unit covariance (O(n^2) memory). Use variance = 'boot' for point-GP inference.",
            call. = FALSE)
  if (variance == "sandwich") obj$vcov <- sandwich_vcov(obj)
  else if (variance == "boot") obj$vcov <- bootstrap_vcov(obj, data)
  obj
}

#' Resolve the CAR weights / precision / constraint-absorbed incidence.
#'
#' Builds the `spq_weights`, the unit map (observation -> unit), the full incidence,
#' the precision, and the per-component sum-to-zero constraint absorption.
#' @keywords internal
resolve_car <- function(spatial_W, geo, car, car_alpha, n) {
  if (is.null(geo))
    stop("`spatial_error = TRUE` needs `coords` (point) or `areal` (areal).",
         call. = FALSE)

  ## ---- build / validate the weights object first (it defines the unit ordering) ----
  if (inherits(spatial_W, "spq_weights") && identical(spatial_W$type, "nngp") &&
      geo$mode == "areal")
    stop("type='nngp' weights are for point data; use coords (not areal).", call. = FALSE)
  if (inherits(spatial_W, "spq_weights")) {
    spqw <- spatial_W
  } else if (!is.null(spatial_W)) {
    spqw <- spq_weights(spatial_W, type = "supplied")
  } else {
    if (geo$mode == "areal") {
      if (is.null(geo$areal))
        stop("Auto-building CAR weights for areal data needs `areal` (an nb/listw), ",
             "or supply `spatial_W`.", call. = FALSE)
      spqw <- spq_weights(geo$areal, type = "queen")
    } else {
      ucc <- geo$coords[!duplicated(
        apply(geo$coords, 1L, function(z) paste(z, collapse = "_"))), , drop = FALSE]
      spqw <- spq_weights(ucc, type = "knn", k = min(5L, nrow(ucc) - 1L))
    }
  }
  L <- nrow(spqw$W)
  w_ids <- spqw$ids
  if (is.null(w_ids)) w_ids <- as.character(seq_len(L))
  ids <- w_ids

  ## ---- observation -> unit index, aligned to the weights' unit ordering ----
  if (geo$mode == "areal") {
    region <- as.character(geo$region)
    unit_idx <- match(region, w_ids)
    if (anyNA(unit_idx))
      stop("Region labels in `coords` not found among the weights' unit ids. ",
           "Ensure spatial_W$ids match levels(factor(region)).", call. = FALSE)
  } else {
    cc <- geo$coords
    key <- apply(cc, 1L, function(z) paste(z, collapse = "_"))
    ukey <- unique(key)
    if (length(ukey) != L)
      stop(sprintf("point data has %d distinct locations but spatial_W has %d units.",
                   length(ukey), L), call. = FALSE)
    unit_idx <- match(key, ukey)
    ## per-unit training coordinates, aligned to the unit ordering (ukey order), so
    ## predict() can match new point coords to the nearest training unit.
    unit_coords <- cc[match(ukey, key), , drop = FALSE]
    rownames(unit_coords) <- NULL
  }

  is_nngp <- identical(spqw$type, "nngp")
  ## NNGP: proper full-rank precision -> use Q directly, single 'component', NO
  ## sum-to-zero (the proper precision + intercept identify phi; Datta et al. 2016).
  Q <- if (is_nngp) spqw$Q else make_car_precision(spqw, alpha = car_alpha, car = car)
  R <- incidence_matrix(unit_idx, L)
  membership <- if (is_nngp) rep(1L, L) else components_from_W(spqw$W)$membership
  ab <- absorb_car_constraint(R, Q, membership, constrain = !is_nngp)
  list(spqw = spqw, Q = Q, R = R, Rt = ab$Rt, Qt = ab$Qt, Tmat = ab$Tmat,
       Lp = ab$Lp, unit_idx = unit_idx, ids = ids, membership = membership,
       mode = geo$mode, kind = if (is_nngp) "nngp" else "car",
       unit_coords = if (geo$mode == "point") unit_coords else NULL)
}

#' Assemble the obj$car slot from the fitted reduced CAR coefficients.
#' @keywords internal
build_car_slot <- function(beta, car_block, car_obj, car, car_alpha, lam_phi) {
  G <- ncol(beta); L <- nrow(car_obj$Q)
  phi <- matrix(0, L, G)
  if (length(car_block) > 0L) {
    red <- beta[car_block, , drop = FALSE]              # L' x G
    phi <- as.matrix(car_obj$Tmat %*% red)              # L x G (mean-zero per comp)
  }
  rownames(phi) <- car_obj$ids
  colnames(phi) <- paste0("regime", seq_len(G))
  list(phi = phi, W = car_obj$spqw, Q = car_obj$Q, alpha = car_alpha, car = car,
       kind = if (!is.null(car_obj$kind)) car_obj$kind else "car",
       lambda = lam_phi, Tmat = car_obj$Tmat, car_block = car_block,
       units = list(ids = car_obj$ids, unit_idx = car_obj$unit_idx,
                    membership = car_obj$membership, mode = car_obj$mode,
                    coords = car_obj$unit_coords))
}

#' Resolve coordinates / areal specification.
#'
#' For the CAR spatial-error path with a supplied `spatial_W` (areal data without an
#' \pkg{spdep} `nb`), `coords` is treated as a length-n vector of region labels even
#' when `areal` is `NULL`.
#' @keywords internal
resolve_coords <- function(coords, areal, data, n, spatial_error = FALSE,
                           spatial_W = NULL) {
  if (!is.null(areal)) {
    region <- if (is.character(coords) && length(coords) == 1L) data[[coords]] else coords
    if (is.null(region)) stop("Areal fit needs `coords` = region labels (length n).")
    return(list(mode = "areal", region = as.factor(region), areal = areal,
                coords = NULL))
  }
  ## CAR-error with a supplied W and region labels (no nb / no mgcv basis): areal mode.
  ## Disambiguate by length, NOT type: a length-2 character vector is the documented
  ## point-coords column-name form (coords = c("sx","sy")), not two region labels; and a
  ## 2-column matrix is always point coords. Region labels are either a single column
  ## name (length 1) or a length-n vector/factor.
  is_colname_pair <- is.character(coords) && length(coords) == 2L &&
    all(coords %in% names(data))
  is_point_matrix <- (is.matrix(coords) || is.data.frame(coords)) && ncol(coords) == 2L
  if (isTRUE(spatial_error) && !is.null(spatial_W) && !is.null(coords) &&
      !is_colname_pair && !is_point_matrix &&
      (is.character(coords) || is.factor(coords) ||
       (is.atomic(coords) && !is.matrix(coords) && length(coords) == n))) {
    region <- if (is.character(coords) && length(coords) == 1L) data[[coords]] else coords
    return(list(mode = "areal", region = as.factor(region), areal = NULL,
                coords = NULL))
  }
  if (is.null(coords)) stop("Spatial terms need `coords` (point) or `areal` (areal).")
  cc <- if (is.character(coords)) as.matrix(data[, coords, drop = FALSE]) else as.matrix(coords)
  if (ncol(cc) != 2L) stop("`coords` must have two columns for point data.")
  storage.mode(cc) <- "double"
  list(mode = "point", coords = cc, areal = NULL, region = NULL)
}

#' Assemble augmented component / gate designs and their penalties.
#'
#' When `car` is non-`NULL` (the CAR spatial-error path) it must be a list with the
#' constraint-absorbed incidence `Rt` (`n x L'`), reduced precision `Qt` (`L' x L'`),
#' and `lambda` (the CAR penalty weight `lambda_phi`). The CAR columns are appended to
#' the component design `Xt` and the penalty `lambda * Qt` is placed on the CAR block;
#' `car_block` records the CAR column indices (so relabelling / ordering ignore them).
#' @keywords internal
build_designs <- function(X, W, B, Omega, slope_idx, spatial_gate, spatial_coef,
                          lam_g, lam_b, gate_ridge, r, car = NULL) {
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
  ## ---- append the constraint-absorbed CAR block (spatial-error path) ----
  car_block <- integer(0); car_sparse <- FALSE
  if (!is.null(car)) {
    Lp <- car$Lp
    if (Lp > 0L) {
      Pbeta_full <- Pen_beta
      ## assemble sparse augmented design and penalty
      Xt_s <- methods::as(methods::as(Matrix::Matrix(Xt, sparse = TRUE),
                                      "CsparseMatrix"), "generalMatrix")
      Xt <- cbind(Xt_s, car$Rt)
      car_block <- (P + 1L):(P + Lp)
      Ptot <- P + Lp
      Pen_beta <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                       x = numeric(0), dims = c(Ptot, Ptot))
      if (any(Pbeta_full != 0)) {
        nz <- which(Pbeta_full != 0, arr.ind = TRUE)
        Pen_beta <- Pen_beta +
          Matrix::sparseMatrix(i = nz[, 1L], j = nz[, 2L],
                               x = Pbeta_full[nz], dims = c(Ptot, Ptot))
      }
      Pen_beta[car_block, car_block] <- car$lambda * car$Qt
      Pen_beta <- methods::as(Pen_beta, "CsparseMatrix")
      car_sparse <- TRUE
    }
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
       gate_spat = gate_spat, qw = ncol(W), r = r, spatial_coef = spatial_coef,
       car_block = car_block, car_sparse = car_sparse)
}

#' Constant coefficient matrix (p x G) from the augmented coefficients.
#'
#' Returns only the beta (intercept + constant-slope) rows; the spatial-slope basis
#' columns and the CAR block (`des$car_block`) are excluded, so relabelling, ordering,
#' and label alignment never see the high-dimensional CAR noise.
#' @keywords internal
extract_const <- function(beta, des, p) {
  cr <- des$const_rows
  if (length(des$car_block %||% integer(0)) == 0L &&
      length(cr) == nrow(beta)) return(beta)
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
X_fitted <- function(Xt, beta) as.matrix(Xt %*% beta)

#' Effective degrees of freedom (component smoothers + gate smoother).
#'
#' The component trace `tr[(X'WpX + P)^{-1} X'WpX]` is computed sparsely when the CAR
#' block is present (the design is a sparse `Matrix`), adding the CAR-smoother trace.
#' @keywords internal
compute_edf <- function(fit, Xt, Z, Pen_beta, Pen_gamma, G, des, spatial_coef, w_obs = NULL) {
  comp <- 0
  if (is.null(w_obs)) w_obs <- rep(1, nrow(Xt))
  car_on <- isTRUE(des$car_sparse)
  for (k in seq_len(G)) {
    e <- w_obs * fit$posterior[, k]                  # weight ~ obs weight * responsibility
    if (car_on) {
      A <- Matrix::crossprod(Xt, e * Xt)
      H <- A + Pen_beta + 1e-8 * Matrix::Diagonal(ncol(Xt))
      comp <- comp + tryCatch(sum(Matrix::diag(Matrix::solve(H, A))),
                              error = function(er) ncol(Xt))
    } else {
      A <- crossprod(Xt, e * Xt)
      H <- A + Pen_beta + 1e-8 * diag(ncol(Xt))
      comp <- comp + sum(diag(safe_solve(H, A)))
    }
  }
  ## spatial portion = component edf beyond the constant (intercept + flat-slope) part
  has_surface <- length(des$spat_blocks) > 0L || length(des$car_block %||% integer(0)) > 0L
  n_const <- if (has_surface) G * length(des$const_rows) else 0
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

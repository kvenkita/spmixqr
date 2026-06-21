#' Select regimes and smoothing parameters
#'
#' Fits [spmixqr] over a grid of `G` and/or roughness penalties and returns the best
#' by BIC (post-hoc unsmoothed log-likelihood plus effective df) or by K-fold
#' held-out predictive check loss. For the KDE density path the log-likelihood is
#' unavailable, so `criterion = "cv"` (the check-loss surrogate) is used.
#'
#' @param formula,data,coords,areal,tau,gating,spatial_gate,spatial_coef,method
#'   Passed to [spmixqr()].
#' @param spatial_error,spatial_W,car,car_alpha CAR spatial-error settings passed to
#'   [spmixqr()]. When `spatial_error = TRUE` the `lambda_error_grid` is searched
#'   (BIC, with the CAR-effective-df term) and `spatial_gate` is forced off (the
#'   guardrail); `car_alpha` is fixed (not selected by the check loss).
#' @param G_grid Candidate regime counts.
#' @param lambda_gate_grid,lambda_coef_grid,lambda_error_grid Candidate penalties.
#' @param criterion `"bic"` or `"cv"`.
#' @param folds Number of CV folds (for `criterion = "cv"`).
#' @param control A [spmixqr_control()] list.
#' @return A list with the best fit, the chosen settings, and the score table.
#' @export
spmixqr_select <- function(formula, data, coords = NULL, areal = NULL, tau = 0.5,
                           gating = ~1, spatial_gate = TRUE, spatial_coef = TRUE,
                           spatial_error = FALSE, spatial_W = NULL,
                           car = c("proper", "icar"), car_alpha = 0.95,
                           method = c("ald", "kde"), G_grid = 2:3,
                           lambda_gate_grid = c(0.1, 1, 10),
                           lambda_coef_grid = c(0.1, 1, 10),
                           lambda_error_grid = c(0.1, 1, 10),
                           criterion = c("bic", "cv"), folds = 5L,
                           control = spmixqr_control()) {
  method <- match.arg(method); criterion <- match.arg(criterion); car <- match.arg(car)
  if (method == "kde") criterion <- "cv"
  ## under the guardrail spatial_error forces spatial_gate off (pass FALSE explicitly
  ## so the internal re-calls never trip the guardrail or re-message B times).
  sg <- if (isTRUE(spatial_error)) FALSE else spatial_gate
  le_grid <- if (isTRUE(spatial_error)) lambda_error_grid else NA
  grid <- expand.grid(G = G_grid, lambda_gate = lambda_gate_grid,
                      lambda_coef = lambda_coef_grid, lambda_error = le_grid,
                      KEEP.OUT.ATTRS = FALSE)
  n <- nrow(data)
  fold_id <- if (criterion == "cv") sample(rep_len(seq_len(folds), n)) else NULL

  fit_one <- function(dat, cc, g, lg, lc, le) {
    spmixqr(formula, dat, cc, areal, G = g, tau = tau, gating = gating,
            spatial_gate = sg, spatial_coef = spatial_coef,
            spatial_error = spatial_error, spatial_W = spatial_W,
            car = car, car_alpha = car_alpha,
            method = method, lambda_gate = lg, lambda_coef = lc,
            lambda_error = if (isTRUE(spatial_error)) le else NULL,
            variance = "none", control = control)
  }

  score_row <- function(g, lg, lc, le) {
    if (criterion == "bic") {
      fit <- fit_one(data, coords, g, lg, lc, le)
      return(fit$bic)
    }
    err <- 0
    for (f in seq_len(folds)) {
      tr <- fold_id != f; te <- !tr
      cc_tr <- coords_subset(coords, tr); cc_te <- coords_subset(coords, te)
      fit <- tryCatch(fit_one(data[tr, , drop = FALSE], cc_tr, g, lg, lc, le),
                      error = function(e) NULL)
      if (is.null(fit)) { err <- err + Inf; next }
      qhat <- tryCatch(predict(fit, newdata = data[te, , drop = FALSE],
                               newcoords = cc_te, type = "quantile"),
                       error = function(e) NULL)
      if (is.null(qhat)) { err <- err + Inf; next }
      yte <- data[te, all.vars(formula)[1L]]
      err <- err + sum(rho_tau(yte - qhat, tau))
    }
    err
  }

  grid$score <- mapply(score_row, grid$G, grid$lambda_gate, grid$lambda_coef,
                       grid$lambda_error)
  best <- grid[which.min(grid$score), , drop = FALSE]
  fit <- fit_one(data, coords, best$G, best$lambda_gate, best$lambda_coef,
                 best$lambda_error)
  list(fit = fit, best = best, criterion = criterion, table = grid)
}

#' Subset coordinates (matrix or region vector) by a logical index.
#' @keywords internal
coords_subset <- function(coords, idx) {
  if (is.null(coords)) return(NULL)
  if (is.matrix(coords) || is.data.frame(coords)) coords[idx, , drop = FALSE]
  else coords[idx]
}

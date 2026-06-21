#' Select regimes and smoothing parameters
#'
#' Fits [spmixqr] over a grid of `G` and/or roughness penalties and returns the best
#' by BIC (post-hoc unsmoothed log-likelihood plus effective df) or by K-fold
#' held-out predictive check loss. For the KDE density path the log-likelihood is
#' unavailable, so `criterion = "cv"` (the check-loss surrogate) is used.
#'
#' @param formula,data,coords,areal,tau,gating,spatial_gate,spatial_coef,method
#'   Passed to [spmixqr()].
#' @param G_grid Candidate regime counts.
#' @param lambda_gate_grid,lambda_coef_grid Candidate penalties.
#' @param criterion `"bic"` or `"cv"`.
#' @param folds Number of CV folds (for `criterion = "cv"`).
#' @param control A [spmixqr_control()] list.
#' @return A list with the best fit, the chosen settings, and the score table.
#' @export
spmixqr_select <- function(formula, data, coords = NULL, areal = NULL, tau = 0.5,
                           gating = ~1, spatial_gate = TRUE, spatial_coef = TRUE,
                           method = c("ald", "kde"), G_grid = 2:3,
                           lambda_gate_grid = c(0.1, 1, 10),
                           lambda_coef_grid = c(0.1, 1, 10),
                           criterion = c("bic", "cv"), folds = 5L,
                           control = spmixqr_control()) {
  method <- match.arg(method); criterion <- match.arg(criterion)
  if (method == "kde") criterion <- "cv"
  grid <- expand.grid(G = G_grid, lambda_gate = lambda_gate_grid,
                      lambda_coef = lambda_coef_grid, KEEP.OUT.ATTRS = FALSE)
  n <- nrow(data)
  fold_id <- if (criterion == "cv") sample(rep_len(seq_len(folds), n)) else NULL

  score_row <- function(g, lg, lc) {
    if (criterion == "bic") {
      fit <- spmixqr(formula, data, coords, areal, G = g, tau = tau, gating = gating,
                     spatial_gate = spatial_gate, spatial_coef = spatial_coef,
                     method = method, lambda_gate = lg, lambda_coef = lc,
                     variance = "none", control = control)
      return(fit$bic)
    }
    err <- 0
    for (f in seq_len(folds)) {
      tr <- fold_id != f; te <- !tr
      cc_tr <- coords_subset(coords, tr); cc_te <- coords_subset(coords, te)
      fit <- tryCatch(spmixqr(formula, data[tr, , drop = FALSE], cc_tr, areal,
                              G = g, tau = tau, gating = gating,
                              spatial_gate = spatial_gate, spatial_coef = spatial_coef,
                              method = method, lambda_gate = lg, lambda_coef = lc,
                              variance = "none", control = control),
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

  grid$score <- mapply(score_row, grid$G, grid$lambda_gate, grid$lambda_coef)
  best <- grid[which.min(grid$score), , drop = FALSE]
  fit <- spmixqr(formula, data, coords, areal, G = best$G, tau = tau, gating = gating,
                 spatial_gate = spatial_gate, spatial_coef = spatial_coef,
                 method = method, lambda_gate = best$lambda_gate,
                 lambda_coef = best$lambda_coef, control = control)
  list(fit = fit, best = best, criterion = criterion, table = grid)
}

#' Subset coordinates (matrix or region vector) by a logical index.
#' @keywords internal
coords_subset <- function(coords, idx) {
  if (is.null(coords)) return(NULL)
  if (is.matrix(coords) || is.data.frame(coords)) coords[idx, , drop = FALSE]
  else coords[idx]
}

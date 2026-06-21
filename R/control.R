#' Control parameters for [spmixqr()]
#'
#' Tuning knobs for the penalised EM, the spatial basis, label handling, and the
#' bootstrap. All have sensible defaults.
#'
#' @param nstart Number of EM starts (the mixture likelihood is multimodal; Wu &
#'   Yao 2016 stress initial-value sensitivity). The best start by penalised
#'   smoothed objective is kept.
#' @param maxit,tol EM iteration cap and relative-change tolerance.
#' @param basis_type Spatial basis: `"tp"` (thin-plate), `"gp"` (low-rank Gaussian
#'   process) for point data, or `"mrf"` (Markov random field) for areal data.
#' @param k Spatial-basis dimension (clamped to the number of unique
#'   locations/regions minus one).
#' @param scale_coords Standardise coordinates before building the basis (kept as a
#'   stored transform for prediction at new locations). Recommended `TRUE`.
#' @param label_order Component ordering key for relabelling: `"slope"` (default,
#'   per-regime average ordering-covariate slope) or `"intercept"`.
#' @param order_var Index (among the slope covariates) used as the ordering
#'   covariate; default `1` (the first non-intercept covariate).
#' @param gate_ridge Ridge on the non-spatial gate coefficients (also stabilises the
#'   ICAR/MRF null space). Matches the `mixqrgate` default for the no-spatial reduction.
#' @param gate_maxit,gate_tol Inner gate Newton/IRLS controls.
#' @param sm_scale,sm_floor Smoothing-bandwidth tuning constant and floor for the
#'   convolution-smoothed component M-step.
#' @param coef_maxit,coef_tol Inner penalised-Newton controls for the component step.
#' @param bandwidth,kde_grid Passed to [mixqr::mixqr_control()] for the KDE density path.
#' @param boot_B,boot_block Bootstrap replicates and number of spatial blocks (per
#'   axis) for the spatial-block bootstrap.
#' @param min_sites Minimum distinct locations required before a block bootstrap is
#'   attempted (it under-covers with few blocks; Lahiri 2003).
#' @param trace Logical; print EM progress.
#' @param seed Optional RNG seed (honoured throughout for reproducibility).
#' @return A list of control parameters.
#' @export
spmixqr_control <- function(nstart = 5L, maxit = 200L, tol = 1e-5,
                            basis_type = c("tp", "gp", "mrf"), k = 20L,
                            scale_coords = TRUE,
                            label_order = c("slope", "intercept"), order_var = 1L,
                            gate_ridge = 1e-3, gate_maxit = 50L, gate_tol = 1e-8,
                            sm_scale = 1, sm_floor = 1e-3,
                            coef_maxit = 50L, coef_tol = 1e-7,
                            bandwidth = NULL, kde_grid = 512L,
                            boot_B = 200L, boot_block = 4L, min_sites = 12L,
                            trace = FALSE, seed = NULL) {
  basis_type <- match.arg(basis_type)
  label_order <- match.arg(label_order)
  list(nstart = as.integer(nstart), maxit = as.integer(maxit), tol = tol,
       basis_type = basis_type, k = as.integer(k), scale_coords = scale_coords,
       label_order = label_order, order_var = as.integer(order_var),
       gate_ridge = gate_ridge, gate_maxit = as.integer(gate_maxit),
       gate_tol = gate_tol, sm_scale = sm_scale, sm_floor = sm_floor,
       coef_maxit = as.integer(coef_maxit), coef_tol = coef_tol,
       bandwidth = bandwidth, kde_grid = as.integer(kde_grid),
       boot_B = as.integer(boot_B), boot_block = as.integer(boot_block),
       min_sites = as.integer(min_sites), trace = isTRUE(trace), seed = seed)
}

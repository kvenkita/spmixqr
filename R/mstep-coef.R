## Penalised spatially-varying-slope component M-step. Minimises the
## convolution-smoothed, ridge/roughness-penalised weighted check loss
##   f(b) = sum_i w_i rho_{tau,h}(y_i - Xtilde_i' b) + (1/2) b' Pen b
## by penalised Newton. As h -> 0 (and Pen = 0) this recovers exact weighted
## quantile regression (mixqr::weighted_rq). Penalty convention matches eq. (4):
## Pen carries lambda_coef*Omega on each spatial-slope block, 0 on the scalar
## intercept and constant-slope columns.

#' Penalised smoothed weighted quantile regression (one regime's component step).
#'
#' Solves the penalised smoothed objective by penalised Newton with **bandwidth
#' annealing** (a wide bandwidth is reduced geometrically to the target, warm-
#' starting each stage) and step-halving. Annealing is essential: at a tiny target
#' bandwidth the kernel Hessian is near-singular and undamped Newton diverges
#' (verified). At convergence the solution approximates the target-bandwidth
#' minimiser, and as the target bandwidth shrinks it approaches exact weighted
#' quantile regression.
#'
#' @param Xt Augmented component design `n x P` (`[1, x_j, x_j*B(s)]`, j>=2).
#' @param y Response.
#' @param tau Quantile level.
#' @param w Observation weights (responsibilities) length n.
#' @param Pen `P x P` penalty matrix (PSD).
#' @param h Bandwidth *rate* (dimensionless); the absolute bandwidth is this rate
#'   times the residual scale, floored at `floor_abs`.
#' @param floor_abs Absolute minimum bandwidth.
#' @param beta_init Optional warm start (ignored if `NULL` or all-zero).
#' @param maxit,tol Newton iterations / tolerance per annealing stage.
#' @return list(beta, hessian, fitted, h (absolute bandwidth used), converged).
#' @keywords internal
pen_smooth_wqr <- function(Xt, y, tau, w, Pen, h, floor_abs = 1e-3,
                           beta_init = NULL, maxit = 50L, tol = 1e-7) {
  P <- ncol(Xt)
  ## warm start: penalised weighted LS (full rank thanks to Pen + tiny ridge).
  ## Always used when no real warm start is supplied (an all-zero start would put
  ## residuals at the response level, collapsing the kernel Hessian).
  b <- if (is.null(beta_init) || all(beta_init == 0)) {
    A <- crossprod(Xt, w * Xt) + Pen + 1e-6 * diag(P)
    as.numeric(safe_solve(A, crossprod(Xt, w * y)))
  } else beta_init
  ## scale-equivariant target bandwidth from the warm-start residual spread
  e0 <- as.numeric(y - Xt %*% b)
  s0 <- stats::mad(e0); if (!is.finite(s0) || s0 <= 0) s0 <- stats::sd(e0)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  h_eff <- max(floor_abs, h * s0)
  ## anneal from a bandwidth wide enough to cover the residuals down to h_eff
  h_start <- max(h_eff, 1.5 * s0)
  nstage <- max(1L, ceiling(log(h_start / h_eff) / log(2)) + 1L)
  h_path <- if (nstage == 1L) h_eff
            else exp(seq(log(h_start), log(h_eff), length.out = nstage))
  h <- h_eff
  H <- NULL
  for (hk in h_path) {
    for (it in seq_len(maxit)) {
      e <- as.numeric(y - Xt %*% b)
      psi <- psi_smooth(e, tau, hk)
      kw <- w * k_smooth(e, hk)
      grad <- crossprod(Xt, w * psi) - Pen %*% b
      H <- crossprod(Xt, kw * Xt) + Pen + 1e-8 * diag(P)
      step <- as.numeric(safe_solve(H, grad))
      ## step-halving on the smoothed objective
      f0 <- coef_objective(b, Xt, y, tau, w, Pen, hk)
      s <- 1; ok <- FALSE
      for (ls in 1:20) {
        bn <- b + s * step
        if (coef_objective(bn, Xt, y, tau, w, Pen, hk) <= f0 + 1e-12) { ok <- TRUE; break }
        s <- s / 2
      }
      if (!ok) bn <- b + step
      done <- max(abs(bn - b)) < tol
      b <- bn
      if (done) break
    }
  }
  list(beta = as.numeric(b), hessian = H,
       fitted = as.numeric(Xt %*% b), h = h, converged = TRUE)
}

#' Sparse penalised smoothed weighted quantile regression for the CAR design.
#'
#' Identical annealed penalised-Newton logic to [pen_smooth_wqr()], but the augmented
#' design `[X, Rt]` (with `Rt` the constraint-absorbed CAR incidence) and the penalty
#' block are stored as \pkg{Matrix} sparse objects, and the `(p + L')^2` Newton system
#' is assembled and solved sparsely via [Matrix::solve()] / a sparse Cholesky. The
#' dense [pen_smooth_wqr()] stays the non-CAR path; this solver is the CAR path
#' (`spatial_error = TRUE`). A dense fallback is used for tiny `L'`.
#'
#' @param Xt Augmented component design, a sparse `Matrix` `n x P` (`P = p + L'`).
#' @param y Response.
#' @param tau Quantile level.
#' @param w Observation weights (responsibilities) length n.
#' @param Pen `P x P` sparse penalty matrix (PSD; zero on the beta block).
#' @param h Bandwidth rate.
#' @param floor_abs Absolute minimum bandwidth.
#' @param beta_init Optional warm start.
#' @param maxit,tol Newton iterations / tolerance per annealing stage.
#' @param dense_max If `P <= dense_max`, fall back to the dense solver path.
#' @return list(beta, hessian, fitted, h, converged) (beta a plain numeric vector).
#' @keywords internal
pen_smooth_wqr_car <- function(Xt, y, tau, w, Pen, h, floor_abs = 1e-3,
                               beta_init = NULL, maxit = 50L, tol = 1e-7,
                               dense_max = 40L) {
  P <- ncol(Xt)
  if (P <= dense_max)
    return(pen_smooth_wqr(as.matrix(Xt), y, tau, w, as.matrix(Pen), h,
                          floor_abs = floor_abs, beta_init = beta_init,
                          maxit = maxit, tol = tol))
  Xt <- methods::as(methods::as(Xt, "CsparseMatrix"), "generalMatrix")
  Pen <- methods::as(Pen, "CsparseMatrix")
  ridgeP <- function(eps) eps * Matrix::Diagonal(P)
  sp_solve <- function(A, b) {
    out <- tryCatch(as.numeric(Matrix::solve(A, b)), error = function(e) NULL)
    if (is.null(out))
      out <- as.numeric(Matrix::solve(A + ridgeP(1e-6), b))
    out
  }
  ## warm start: penalised weighted LS (sparse)
  b <- if (is.null(beta_init) || all(beta_init == 0)) {
    A <- Matrix::crossprod(Xt, w * Xt) + Pen + ridgeP(1e-6)
    sp_solve(A, Matrix::crossprod(Xt, w * y))
  } else as.numeric(beta_init)
  e0 <- as.numeric(y - Xt %*% b)
  s0 <- stats::mad(e0); if (!is.finite(s0) || s0 <= 0) s0 <- stats::sd(e0)
  if (!is.finite(s0) || s0 <= 0) s0 <- 1
  h_eff <- max(floor_abs, h * s0)
  h_start <- max(h_eff, 1.5 * s0)
  nstage <- max(1L, ceiling(log(h_start / h_eff) / log(2)) + 1L)
  h_path <- if (nstage == 1L) h_eff
            else exp(seq(log(h_start), log(h_eff), length.out = nstage))
  h <- h_eff
  H <- NULL
  for (hk in h_path) {
    for (it in seq_len(maxit)) {
      e <- as.numeric(y - Xt %*% b)
      psi <- psi_smooth(e, tau, hk)
      kw <- w * k_smooth(e, hk)
      grad <- as.numeric(Matrix::crossprod(Xt, w * psi) - Pen %*% b)
      H <- Matrix::crossprod(Xt, kw * Xt) + Pen + ridgeP(1e-8)
      step <- sp_solve(H, grad)
      f0 <- coef_objective_sp(b, Xt, y, tau, w, Pen, hk)
      s <- 1; ok <- FALSE
      for (ls in 1:20) {
        bn <- b + s * step
        if (coef_objective_sp(bn, Xt, y, tau, w, Pen, hk) <= f0 + 1e-12) { ok <- TRUE; break }
        s <- s / 2
      }
      if (!ok) bn <- b + step
      done <- max(abs(bn - b)) < tol
      b <- bn
      if (done) break
    }
  }
  list(beta = as.numeric(b), hessian = H,
       fitted = as.numeric(Xt %*% b), h = h, converged = TRUE)
}

#' Sparse-aware smoothed penalised objective (Matrix-friendly).
#' @keywords internal
coef_objective_sp <- function(b, Xt, y, tau, w, Pen, h) {
  e <- as.numeric(y - Xt %*% b)
  sum(w * rho_smooth(e, tau, h)) + 0.5 * as.numeric(Matrix::crossprod(b, Pen %*% b))
}

#' Smoothed penalised objective value for one regime (for tests and the EM objective).
#' @keywords internal
coef_objective <- function(b, Xt, y, tau, w, Pen, h) {
  e <- as.numeric(y - Xt %*% b)
  sum(w * rho_smooth(e, tau, h)) + 0.5 * as.numeric(crossprod(b, Pen %*% b))
}

#' Classification-conditional smoothed-QR sandwich for one regime's coefficients.
#' Bread-meat-bread form with the density read off the fitted smoothed kernel (never
#' an ALD stand-in). Disclosed as classification-conditional.
#'
#' On the CAR spatial-error path `Xt`/`Pen` arrive as sparse \pkg{Matrix} `dgCMatrix`
#' objects; the bread and meat are dense `(p + L')^2` blocks, so we densify the inputs
#' (`as.matrix`) and work with base linear algebra. This is acceptable because the
#' sandwich already materialises a dense `(p + L')^2` inverse; for very large `L'`
#' (areal L up to ~400 is fine) the bootstrap is the recommended inference path.
#' @keywords internal
coef_sandwich_vcov <- function(Xt, y, tau, w, beta, Pen, h) {
  Xt <- as.matrix(Xt); Pen <- as.matrix(Pen); beta <- as.numeric(beta)
  e <- as.numeric(y - Xt %*% beta)
  psi <- psi_smooth(e, tau, h)
  kw <- w * k_smooth(e, h)
  H <- crossprod(Xt, kw * Xt) + Pen
  meat <- crossprod(Xt, (w^2 * psi^2) * Xt)
  Hinv <- tryCatch(solve(H), error = function(ee) ginv_small(H))
  symmetrise(Hinv %*% meat %*% Hinv)
}

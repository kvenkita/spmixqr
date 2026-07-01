## Spatially penalised multinomial-logit gate. Vendored and generalised from the
## mixqrgate internal `irls_multinom_fit` (Grun & Leisch 2008 / McLachlan & Peel
## 2000 gate Q-function): the scalar ridge `lambda*I` is replaced by a structured
## per-class penalty matrix `Pen` (ridge on the covariate block, `lambda_gate*Omega`
## on the spatial block). Penalty convention matches eq. (4): objective term
## (1/2) sum_k gamma_k' Pen gamma_k, so the gradient carries `- Pen %*% gamma_k` and
## the Hessian `- Pen`.

#' Row-wise softmax over (G-1) linear predictors plus a zero reference column.
#' @keywords internal
softmax_rows <- function(eta) {
  full <- cbind(0, eta)
  full <- full - apply(full, 1L, max)        # overflow guard
  ex <- exp(full)
  ex / rowSums(ex)
}

#' Gate probabilities pi (n x G) from gate coefficients and design.
#' @keywords internal
gate_predict <- function(gamma, Z) softmax_rows(Z %*% gamma)

#' Fit a penalised weighted multinomial logit by Newton/IRLS.
#'
#' @param Z Gate design `n x q1`.
#' @param P Fractional responses `n x G` (responsibilities; rows sum to 1).
#' @param Pen Per-class penalty matrix `q1 x q1` (PSD).
#' @param maxit,tol Newton iterations / tolerance.
#' @return list(gamma `q1 x (G-1)`, hessian, pi, converged).
#' @keywords internal
pen_irls_multinom <- function(Z, P, Pen, maxit = 50L, tol = 1e-8, w = NULL) {
  n <- nrow(Z); q1 <- ncol(Z); G <- ncol(P)
  if (is.null(w)) w <- rep(1, n)
  if (G < 2L)
    return(list(gamma = matrix(0, q1, 0), hessian = matrix(0, 0, 0),
                pi = matrix(1, n, 1), converged = TRUE))
  K <- G - 1L
  gamma <- matrix(0, q1, K)
  pidx <- 2:G
  it <- 0L
  for (it in seq_len(maxit)) {
    pis <- softmax_rows(Z %*% gamma)
    g <- numeric(q1 * K)
    for (a in seq_len(K)) {
      r <- w * (P[, pidx[a]] - pis[, pidx[a]])
      g[((a - 1L) * q1 + 1L):(a * q1)] <- crossprod(Z, r) - Pen %*% gamma[, a]
    }
    H <- matrix(0, q1 * K, q1 * K)
    for (a in seq_len(K)) {
      ia <- ((a - 1L) * q1 + 1L):(a * q1)
      for (b in seq_len(K)) {
        ib <- ((b - 1L) * q1 + 1L):(b * q1)
        wv <- w * pis[, pidx[a]] * ((a == b) - pis[, pidx[b]])
        block <- -crossprod(Z, wv * Z)
        if (a == b) block <- block - Pen
        H[ia, ib] <- block
      }
    }
    step <- safe_solve(H, g)
    gamma_new <- gamma - matrix(step, q1, K)
    if (max(abs(gamma_new - gamma)) < tol) { gamma <- gamma_new; break }
    gamma <- gamma_new
  }
  list(gamma = gamma, hessian = H, pi = softmax_rows(Z %*% gamma),
       converged = it < maxit)
}

#' Penalised gate objective Q(gamma) (for testing the Newton step against optim).
#' @keywords internal
gate_objective <- function(gamma, Z, P, Pen, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(Z))
  pis <- softmax_rows(Z %*% gamma)
  ll <- sum(w * P * log(pmax(pis, .dens_floor)))
  pen <- 0
  for (a in seq_len(ncol(gamma)))
    pen <- pen + 0.5 * as.numeric(crossprod(gamma[, a], Pen %*% gamma[, a]))
  ll - pen
}

#' Classification-conditional sandwich covariance for the gate coefficients.
#' V = Ainv B Ainv, with A the negative penalised-Q Hessian and B the score crossproduct.
#' Independent-score (not spatial-dependence robust) -- disclosed; the bootstrap is
#' the default for spatially-dependent inference.
#' @keywords internal
gate_sandwich_vcov <- function(Z, P, fit, ow = NULL, wtype = "sampling") {
  n <- nrow(Z); q1 <- ncol(Z); G <- ncol(P); K <- G - 1L
  if (K < 1L) return(matrix(0, 0, 0))
  if (is.null(ow)) ow <- rep(1, n)
  pis <- fit$pi; pidx <- 2:G
  S <- matrix(0, n, q1 * K)
  for (a in seq_len(K)) {
    r <- P[, pidx[a]] - pis[, pidx[a]]
    S[, ((a - 1L) * q1 + 1L):(a * q1)] <- Z * r
  }
  wcol <- if (identical(wtype, "sampling")) ow^2 else ow
  A <- -fit$hessian                                # weighted penalised-Q Hessian
  B <- crossprod(S * sqrt(wcol))                   # t(S) diag(wcol) S
  Ainv <- tryCatch(solve(A), error = function(e) ginv_small(A))
  symmetrise(Ainv %*% B %*% Ainv)
}

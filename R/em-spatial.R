## Penalised spatial EM (one start). Reuses mixqr::weighted_rq / constrained_kde for
## the component density machinery and the vendored penalised gate. The component
## M-step is exact weighted QR when slopes are flat, or the smoothed penalised
## solver when slopes vary spatially.

#' One penalised-EM fit from an initial responsibility matrix.
#'
#' When `spatial_error = TRUE` the component M-step is **always** routed through the
#' penalised CAR solver ([pen_smooth_wqr_car()]) on the augmented `[X, Rt]` design with
#' the `lambda_phi * Qt` penalty block, regardless of `spatial_coef` (the unpenalised
#' `weighted_rq` branch would silently drop the CAR penalty for the `G = 1` /
#' constant-slope paths).
#' @keywords internal
spatial_em_fit <- function(y, Xt, Z, G, tau, method, p_init,
                           Pen_beta, Pen_gamma, h, spatial_coef, kdectrl, control,
                           spatial_error = FALSE) {
  n <- length(y); P <- ncol(Xt)
  p <- p_init
  beta <- matrix(0, P, G)
  prev <- NULL; converged <- FALSE; it <- 0L

  h_used <- h
  component_step <- function(p, beta) {
    for (k in seq_len(G)) {
      if (isTRUE(spatial_error)) {
        cf <- pen_smooth_wqr_car(Xt, y, tau, p[, k], Pen_beta, h,
                                 floor_abs = control$sm_floor, beta_init = beta[, k],
                                 maxit = control$coef_maxit, tol = control$coef_tol)
        beta[, k] <- cf$beta
        h_used <<- cf$h
      } else if (spatial_coef) {
        cf <- pen_smooth_wqr(Xt, y, tau, p[, k], Pen_beta, h,
                             floor_abs = control$sm_floor, beta_init = beta[, k],
                             maxit = control$coef_maxit, tol = control$coef_tol)
        beta[, k] <- cf$beta
        h_used <<- cf$h
      } else {
        beta[, k] <- mixqr::weighted_rq(Xt, y, tau, w = p[, k],
          beta_prev = if (any(beta[, k] != 0)) beta[, k] else NULL)
      }
    }
    beta
  }
  densities <- function(beta, p) {
    e <- as.matrix(y - Xt %*% beta)
    if (method == "ald") {
      sigma <- vapply(seq_len(G), function(k) {
        den <- sum(p[, k]); num <- sum(p[, k] * rho_tau(e[, k], tau))
        max(if (den > 0) num / den else 1, 1e-8)
      }, numeric(1))
      dm <- vapply(seq_len(G), function(k) ald_density(e[, k], sigma[k], tau),
                   numeric(n))
      list(sigma = sigma, dens = NULL, dm = dm, e = e)
    } else {
      dens <- mixqr::constrained_kde(e, p, tau, "unequal", kdectrl)
      dm <- vapply(seq_len(G), function(k) pmax(dens[[k]]$eval(e[, k]), .dens_floor),
                   numeric(n))
      list(sigma = NULL, dens = dens, dm = dm, e = e)
    }
  }

  gfit <- NULL
  for (it in seq_len(control$maxit)) {
    beta <- component_step(p, beta)
    cd <- densities(beta, p)
    gfit <- pen_irls_multinom(Z, p, Pen_gamma, control$gate_maxit, control$gate_tol)
    prior <- gfit$pi
    pnew <- normalize_rows(pmax(prior * cd$dm, .dens_floor))
    cur <- c(as.numeric(beta), as.numeric(gfit$gamma))
    if (!is.null(prev)) {
      d <- sum(abs(cur - prev)) / (sum(abs(prev)) + 1e-8)
      if (isTRUE(control$trace)) message(sprintf("  EM iter %3d  rel=%.3e", it, d))
      if (d < control$tol) { converged <- TRUE; p <- pnew; break }
    }
    prev <- cur; p <- pnew
  }

  ## final quantities at the converged responsibilities
  cd <- densities(beta, p)
  if (!spatial_coef) {                       # scale-aware bandwidth for the sandwich
    s0 <- stats::mad(as.numeric(cd$e))
    if (!is.finite(s0) || s0 <= 0) s0 <- stats::sd(as.numeric(cd$e))
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
    h_used <- max(control$sm_floor, h * s0)
  }
  gfit <- pen_irls_multinom(Z, p, Pen_gamma, control$gate_maxit, control$gate_tol)
  prior <- gfit$pi

  pen_b <- 0.5 * sum(vapply(seq_len(G), function(k)
    as.numeric(crossprod(beta[, k], as.numeric(Pen_beta %*% beta[, k]))), numeric(1)))
  pen_g <- if (ncol(gfit$gamma) > 0)
    0.5 * sum(vapply(seq_len(ncol(gfit$gamma)), function(a)
      as.numeric(crossprod(gfit$gamma[, a], Pen_gamma %*% gfit$gamma[, a])), numeric(1)))
    else 0
  loglik <- if (method == "ald")
    sum(log(pmax(rowSums(prior * cd$dm), .dens_floor))) else NA_real_
  objective <- if (is.finite(loglik)) loglik - pen_b - pen_g
               else -sum(p * rho_tau(cd$e, tau)) - pen_b - pen_g

  list(beta = beta, gamma = gfit$gamma, sigma = cd$sigma, dens = cd$dens,
       dm = cd$dm, prior = prior, posterior = p, gate_hessian = gfit$hessian,
       gate_cond = tryCatch(kappa(-gfit$hessian), error = function(e) NA_real_),
       loglik = loglik, objective = objective, pen_beta = pen_b, pen_gamma = pen_g,
       h_used = h_used, iters = it, converged = converged)
}

## S3 method surface, mirroring mixqr / mixqrgate.

#' @export
print.spmixqr <- function(x, ...) {
  cat(sprintf("Spatial mixture of quantile regressions (spmixqr)\n"))
  cat(sprintf("  G = %d regimes,  tau = %.3g,  method = %s\n", x$G, x$tau, x$method))
  cat(sprintf("  spatial gate: %s   spatial slopes: %s   basis: %s (r=%d)\n",
              x$spatial_gate, x$spatial_coef,
              if (!is.null(x$basis)) x$basis$type else "none",
              if (!is.null(x$basis)) x$basis$r else 0L))
  cat("\nConstant component coefficients (intercept + average slopes):\n")
  bc <- x$beta_const
  rownames(bc) <- coef_names(x)
  colnames(bc) <- paste0("regime", seq_len(x$G))
  print(round(bc, 4))
  if (is.finite(x$loglik))
    cat(sprintf("\n  logLik = %.2f   edf = %.2f   AIC = %.1f   BIC = %.1f\n",
                x$loglik, x$edf, x$aic, x$bic))
  if (!is.na(x$diagnostics$label_stability) && x$diagnostics$label_stability > 0.2)
    cat(sprintf("  note: slope surfaces cross at %.0f%% of sites (label coherence approximate).\n",
                100 * x$diagnostics$label_stability))
  invisible(x)
}

#' Coefficient names (intercept + slope variables).
#' @keywords internal
coef_names <- function(x) {
  cn <- colnames(x$X)
  cn[cn == "(Intercept)"] <- "(Intercept)"
  cn
}

#' @export
summary.spmixqr <- function(object, ...) {
  V <- object$vcov
  bc <- object$beta_const
  cn <- coef_names(object)
  ## component constant-coefficient tables (per regime), SEs from the sandwich
  comp_tab <- lapply(seq_len(object$G), function(k) {
    est <- bc[, k]
    se <- rep(NA_real_, length(est))
    if (!is.null(V) && !is.null(V$coef)) {
      cr <- object$design$const_rows
      Vk <- V$coef[[k]]
      if (!is.null(Vk) && all(cr <= nrow(Vk))) se <- sqrt(pmax(diag(Vk)[cr], 0))
    }
    z <- est / se
    data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)), check.names = FALSE,
               row.names = cn)
  })
  names(comp_tab) <- paste0("regime", seq_len(object$G))
  ## gate covariate (alpha) table per non-reference class
  gate_tab <- NULL
  if (object$G > 1L && ncol(object$gamma) > 0L) {
    qw <- object$design$qw; q1 <- nrow(object$gamma)
    wn <- colnames(object$W)
    Vg <- if (!is.null(V) && !is.null(V$gate)) V$gate else NULL
    gate_tab <- lapply(seq_len(ncol(object$gamma)), function(a) {
      est <- object$gamma[seq_len(qw), a]
      se <- rep(NA_real_, qw)
      if (!is.null(Vg)) {
        idx <- (a - 1L) * q1 + seq_len(qw)
        if (all(idx <= nrow(Vg))) se <- sqrt(pmax(diag(Vg)[idx], 0))
      }
      z <- est / se
      data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
                 `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)), check.names = FALSE,
                 row.names = wn)
    })
    names(gate_tab) <- paste0("regime", 2:object$G, " vs regime1")
  }
  ## spatial-error (CAR) block
  car_block <- NULL
  if (isTRUE(object$spatial_error) && !is.null(object$car)) {
    moran_after <- NA_real_; moran_p <- NA_real_
    mo <- object$diagnostics$moran
    if (!is.null(mo)) {
      ms <- if (!is.null(mo$after)) mo$after else mo
      moran_after <- ms$statistic; moran_p <- ms$p_value
    }
    var_share <- apply(object$car$phi, 2L, stats::var)
    kind <- if (!is.null(object$car$kind)) object$car$kind else "car"
    car_block <- list(lambda = object$car$lambda, alpha = object$car$alpha,
                      car = object$car$car, n_units = nrow(object$car$phi),
                      n_comp = object$car$W$n_comp, var_share = var_share,
                      moran = moran_after, moran_p = moran_p, kind = kind,
                      m = object$car$W$m, range = object$car$W$range, nu = object$car$W$nu)
  }
  out <- list(call = object$call, G = object$G, tau = object$tau,
              method = object$method, comp = comp_tab, gate = gate_tab,
              loglik = object$loglik, edf = object$edf, aic = object$aic,
              bic = object$bic, se_method = object$se_method,
              weights_type = if (!is.null(object$prior_weights)) object$weights_type else NULL,
              n_eff = if (!is.null(object$prior_weights)) object$weights_sum else NULL,
              diagnostics = object$diagnostics, spatial_gate = object$spatial_gate,
              spatial_coef = object$spatial_coef, spatial_plus = isTRUE(object$spatial_plus),
              spatial_plus_r2 = object$diagnostics$spatial_plus,
              spatial_error = isTRUE(object$spatial_error), car = car_block)
  class(out) <- "summary.spmixqr"
  out
}

#' @export
print.summary.spmixqr <- function(x, ...) {
  cat("Spatial mixture of quantile regressions (spmixqr)\n")
  cat(sprintf("G = %d,  tau = %.3g,  method = %s\n\n", x$G, x$tau, x$method))
  for (k in seq_along(x$comp)) {
    cat(sprintf("-- Component %s (constant part) --\n", names(x$comp)[k]))
    stats::printCoefmat(x$comp[[k]], P.values = TRUE, has.Pvalue = TRUE, digits = 4)
    cat("\n")
  }
  if (!is.null(x$gate)) {
    cat("-- Spatial gate (covariate log-odds vs regime1) --\n")
    for (a in seq_along(x$gate)) {
      cat(sprintf("   %s:\n", names(x$gate)[a]))
      stats::printCoefmat(x$gate[[a]], P.values = TRUE, has.Pvalue = TRUE, digits = 4)
    }
    cat("\n")
  }
  if (isTRUE(x$spatial_plus) && !is.null(x$spatial_plus_r2)) {
    cat("-- Spatial+ confounding safeguard (Dupont, Wood & Augustin 2022) --\n")
    cat("   covariates residualised against a spatial smooth; slopes below are the effect of the\n")
    cat("   NON-spatial part of each covariate. Mean (least-squares) residualisation deconfounds\n")
    cat("   the slope at EVERY tau under an additive effect (deconfounding is tau-free; see the\n")
    cat("   spatial-error vignette), given the smooth out-resolves the covariate's spatial mean.\n")
    cat("   Spatial R^2 removed:\n")
    r2 <- x$spatial_plus_r2
    for (i in seq_len(nrow(r2)))
      cat(sprintf("     %-16s %.3f\n", r2$variable[i], r2$spatialR2[i]))
    cat("\n")
  }
  if (!is.null(x$car)) {
    cb <- x$car
    if (identical(cb$kind, "nngp")) {
      cat("-- Spatial NNGP error (Datta et al. 2016) --\n")
      cat(sprintf("   units = %d;  Matern nu = %s, range = %.4g (raw units), m = %d neighbours,  lambda_phi = %.4g\n",
                  cb$n_units, as.character(cb$nu), cb$range, cb$m, cb$lambda))
    } else {
      cat("-- Spatial CAR error --\n")
      cat(sprintf("   units = %d (%s connected component(s));  %s CAR, alpha = %.2f (fixed),  lambda_phi = %.4g\n",
                  cb$n_units, as.character(cb$n_comp), cb$car, cb$alpha, cb$lambda))
    }
    cat(sprintf("   per-regime spatial variance: %s\n",
                paste(sprintf("%.3g", cb$var_share), collapse = ", ")))
    if (is.finite(cb$moran))
      cat(sprintf("   residual Moran's I = %.4f (perm p = %.3f)\n", cb$moran, cb$moran_p))
    cat("   note: lambda_phi is selectable by BIC via spmixqr_select() (this fit holds it fixed);\n")
    if (identical(cb$kind, "nngp"))
      cat("         range/nu are weakly identified (fixed/selected, not estimated; Zhang 2004).\n\n")
    else
      cat("         alpha is weakly identified (fixed, not estimated).\n\n")
  }
  cat(sprintf("SE method: %s", x$se_method))
  if (identical(x$se_method, "sandwich"))
    cat(" (classification-conditional; use variance='boot' for reporting)")
  cat("\n")
  if (!is.null(x$weights_type))
    cat(sprintf("weights: %s (effective n = %.1f)\n", x$weights_type, x$n_eff))
  if (is.finite(x$loglik))
    cat(sprintf("logLik %.2f   edf %.2f   AIC %.1f   BIC %.1f\n",
                x$loglik, x$edf, x$aic, x$bic))
  d <- x$diagnostics
  cat(sprintf("converged %s | starts %d | gate cond %.1f | mean class entropy %.2f\n",
              d$converged, d$n_starts, d$gate_cond, d$class_entropy))
  if (x$spatial_coef && !is.na(d$label_stability))
    cat(sprintf("label stability: slope surfaces cross at %.0f%% of sites\n",
                100 * d$label_stability))
  invisible(x)
}

#' @export
coef.spmixqr <- function(object, ...) {
  out <- object$beta_const
  rownames(out) <- coef_names(object)
  colnames(out) <- paste0("regime", seq_len(object$G))
  out
}

#' @export
vcov.spmixqr <- function(object, ...) {
  if (is.null(object$vcov)) object$vcov <- sandwich_vcov(object)
  object$vcov
}

#' @export
confint.spmixqr <- function(object, parm, level = 0.95, ...) {
  V <- vcov.spmixqr(object)
  cn <- coef_names(object); cr <- object$design$const_rows
  za <- stats::qnorm(1 - (1 - level) / 2)
  res <- lapply(seq_len(object$G), function(k) {
    est <- object$beta_const[, k]
    se <- if (!is.null(V$coef) && !is.null(V$coef[[k]]))
      sqrt(pmax(diag(V$coef[[k]])[cr], 0)) else rep(NA_real_, length(est))
    m <- cbind(est - za * se, est + za * se)
    colnames(m) <- c(sprintf("%.1f%%", 100 * (1 - level) / 2),
                     sprintf("%.1f%%", 100 * (1 + level) / 2))
    rownames(m) <- cn; m
  })
  names(res) <- paste0("regime", seq_len(object$G))
  res
}

#' @export
logLik.spmixqr <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$edf
  attr(val, "nobs") <- if (isTRUE(object$weights_type == "frequency") &&
                           !is.null(object$prior_weights)) object$weights_sum
                       else length(object$y)
  class(val) <- "logLik"
  val
}

#' @export
AIC.spmixqr <- function(object, ..., k = 2) {
  if (!is.finite(object$loglik)) return(NA_real_)
  -2 * object$loglik + k * object$edf
}

#' @export
BIC.spmixqr <- function(object, ...) object$bic

#' @export
nobs.spmixqr <- function(object, ...) length(object$y)

#' @export
fitted.spmixqr <- function(object, ...) rowSums(object$posterior * object$fitted_q)

#' @export
residuals.spmixqr <- function(object, ...) object$y - fitted.spmixqr(object)

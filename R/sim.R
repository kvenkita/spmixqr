#' Simulate from a spatial quantile mixture model
#'
#' Generates point (or areal) data from a two-regime spatial mixture of quantile
#' regressions: a smooth spatial gate and spatially varying component slopes (à la
#' Reich, Fuentes & Dunson 2011, eq. 17, inside a mixture). Used for the validation
#' study. The conditional `tau`-quantile is exact by construction (the error's
#' `tau`-quantile is subtracted).
#'
#' @param n Number of observations / distinct locations.
#' @param G Number of regimes (2 or 3 supported).
#' @param tau Quantile level.
#' @param gate_slope Strength of the spatial gate (0 = constant gate; negative
#'   control uses 0).
#' @param coef_slope Strength of the spatial slope variation (0 = flat slopes).
#' @param sep Regime separation (intercept gap).
#' @param sd Error scale.
#' @param crossing If `TRUE`, make slope surfaces cross in space (stresses labels).
#' @param error `"normal"` or `"ald"` component error.
#' @param seed Optional seed.
#' @return A list with `data` (y, x, x2), `coords` (n x 2), and `truth`.
#' @examples
#' d <- sim_spmixqr(n = 150, seed = 1)
#' str(d$truth)
#' @export
sim_spmixqr <- function(n = 200L, G = 2L, tau = 0.5, gate_slope = 2.5,
                        coef_slope = 1.5, sep = 3, sd = 0.7, crossing = FALSE,
                        error = c("normal", "ald"), seed = NULL) {
  error <- match.arg(error)
  if (!is.null(seed)) set.seed(seed)
  s1 <- stats::runif(n); s2 <- stats::runif(n)
  coords <- cbind(s1, s2)
  x <- stats::rnorm(n); x2 <- stats::rnorm(n)

  ## spatial gate: logit of regime>1 increases with s1
  cen <- 2 * (s1 - 0.5)
  if (G == 2L) {
    eta <- cbind(0, -0.3 + gate_slope * cen)
  } else {
    eta <- cbind(0, -0.3 + gate_slope * cen, -0.3 + gate_slope * (2 * (s2 - 0.5)))
  }
  pi_mat <- exp(eta) / rowSums(exp(eta))
  z <- apply(pi_mat, 1L, function(pp) sample.int(length(pp), 1L, prob = pp))

  ## spatially varying slope on x; scalar intercepts (separated)
  intercepts <- (seq_len(G) - 1) * sep
  if (crossing) {
    slope_fun <- function(k, s) (if (k == 1) 1 else -1) * (coef_slope * (2 * s - 1)) +
      (if (k == 1) 1 else 2)
  } else {
    slope_fun <- function(k, s) (if (k == 1) 0.5 else -2) + coef_slope * s * (if (k == 1) 1 else -1)
  }
  q <- numeric(n); slope_true <- numeric(n)
  for (i in seq_len(n)) {
    k <- z[i]
    bx <- slope_fun(k, s1[i])
    slope_true[i] <- bx
    q[i] <- intercepts[k] + bx * x[i] + 0.5 * x2[i]
  }
  eps <- if (error == "ald") rald(n, tau, sd) else stats::rnorm(n, 0, sd)
  eps <- eps - stats::quantile(eps, tau)        # make tau-quantile exactly 0
  y <- q + eps

  list(data = data.frame(y = y, x = x, x2 = x2),
       coords = coords,
       truth = list(z = z, pi = pi_mat, q = q, slope_x = slope_true,
                    intercepts = intercepts, tau = tau, G = G,
                    gate_slope = gate_slope, coef_slope = coef_slope))
}

#' Draw asymmetric-Laplace noise (tau-parametrisation).
#' @keywords internal
rald <- function(n, tau, sigma) {
  u <- stats::runif(n)
  ifelse(u < tau, sigma / (1 - tau) * log(u / tau),
         -sigma / tau * log((1 - u) / (1 - tau)))
}

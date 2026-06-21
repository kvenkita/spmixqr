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
#' @param spatial_error If `TRUE`, generate data on a regular `lattice x lattice`
#'   grid and add a per-regime CAR spatial random effect `phi_k` drawn from the
#'   proper-CAR precision built from the lattice's rook contiguity. The returned list
#'   then gains `spatial_W` (the [spq_weights()] object), and `truth$phi` (the L x G
#'   true CAR surface, each column mean-zero) and `truth$unit` (observation -> unit).
#' @param lattice Side length of the square lattice when `spatial_error = TRUE`
#'   (`L = lattice^2` units); `n` observations are assigned to units (with repeats).
#' @param car_alpha,car_rho Proper-CAR alpha used to build the precision and the
#'   marginal SD scaling of the simulated CAR effect (`car_rho` = phi SD).
#' @param seed Optional seed.
#' @return A list with `data` (y, x, x2), `coords` (n x 2), and `truth`. With
#'   `spatial_error = TRUE` also `spatial_W` and a `region` column / `truth$phi`.
#' @examples
#' d <- sim_spmixqr(n = 150, seed = 1)
#' str(d$truth)
#' @export
sim_spmixqr <- function(n = 200L, G = 2L, tau = 0.5, gate_slope = 2.5,
                        coef_slope = 1.5, sep = 3, sd = 0.7, crossing = FALSE,
                        error = c("normal", "ald"), spatial_error = FALSE,
                        lattice = 8L, car_alpha = 0.95, car_rho = 1.5,
                        seed = NULL) {
  error <- match.arg(error)
  if (!is.null(seed)) set.seed(seed)
  if (isTRUE(spatial_error))
    return(sim_spmixqr_car(n, G, tau, coef_slope, sep, sd, error, lattice,
                           car_alpha, car_rho))
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

#' Simulate from a CAR spatial-error quantile mixture on a lattice.
#'
#' Lattice (`L = lattice^2` units, rook contiguity); each regime carries a proper-CAR
#' spatial random effect drawn from `N(0, [Q(alpha)]^{-1})`, scaled to SD `car_rho`
#' and centred (mean-zero), plus a constant slope on `x`. The conditional tau-quantile
#' is exact by construction. Returns `spatial_W`, `region`, `truth$phi`, `truth$unit`.
#' @keywords internal
sim_spmixqr_car <- function(n, G, tau, coef_slope, sep, sd, error, lattice,
                            car_alpha, car_rho) {
  side <- as.integer(lattice); L <- side * side
  ## rook contiguity weights on the lattice
  coordsL <- expand.grid(r = seq_len(side), c = seq_len(side))
  W <- matrix(0, L, L)
  for (a in seq_len(L)) for (b in seq_len(L)) {
    if (a < b) {
      man <- abs(coordsL$r[a] - coordsL$r[b]) + abs(coordsL$c[a] - coordsL$c[b])
      if (man == 1) { W[a, b] <- 1; W[b, a] <- 1 }
    }
  }
  spqw <- spq_weights(W, type = "supplied", ids = as.character(seq_len(L)))
  Q <- as.matrix(make_car_precision(spqw, alpha = car_alpha, car = "proper"))
  Sig <- solve(Q)
  ## draw per-regime CAR surfaces (mean-zero, scaled to car_rho)
  ev <- eigen(Sig, symmetric = TRUE)
  Aroot <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0))) %*% t(ev$vectors)
  phi <- matrix(0, L, G)
  for (k in seq_len(G)) {
    z <- as.numeric(Aroot %*% stats::rnorm(L))
    z <- z - mean(z)
    z <- z / stats::sd(z) * car_rho
    phi[, k] <- z
  }
  ## assign n observations to units (cycling), draw covariates and regimes
  unit <- sample.int(L, n, replace = TRUE)
  x <- stats::rnorm(n); x2 <- stats::rnorm(n)
  ## constant gate (membership independent of space under the guardrail)
  pi_vec <- rep(1 / G, G)
  z_lat <- sample.int(G, n, replace = TRUE, prob = pi_vec)
  intercepts <- (seq_len(G) - 1) * sep
  slope_x <- if (G == 1L) 0.8 else c(0.8, -1.2, 0.3)[seq_len(G)]
  q <- numeric(n)
  for (i in seq_len(n)) {
    k <- z_lat[i]
    q[i] <- intercepts[k] + slope_x[k] * x[i] + 0.5 * x2[i] + phi[unit[i], k]
  }
  eps <- if (error == "ald") rald(n, tau, sd) else stats::rnorm(n, 0, sd)
  eps <- eps - stats::quantile(eps, tau)
  y <- q + eps
  region <- as.character(unit)
  list(data = data.frame(y = y, x = x, x2 = x2, region = region),
       coords = region, areal = NULL, spatial_W = spqw, region = region,
       truth = list(z = z_lat, phi = phi, unit = unit, q = q,
                    intercepts = intercepts, slope_x = slope_x,
                    tau = tau, G = G, alpha = car_alpha, car_rho = car_rho, L = L))
}

#' Draw asymmetric-Laplace noise (tau-parametrisation).
#' @keywords internal
rald <- function(n, tau, sigma) {
  u <- stats::runif(n)
  ifelse(u < tau, sigma / (1 - tau) * log(u / tau),
         -sigma / tau * log((1 - u) / (1 - tau)))
}

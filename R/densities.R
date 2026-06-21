## Convolution-smoothed check loss (Fernandes, Guerre & Horta 2021; Tan, Wang &
## Zhou 2022, "conquer"). Self-implemented (cited as method, not a dependency).
##
## For residual u and bandwidth h with a Gaussian convolution kernel:
##   rho_{tau,h}(u) = u*(tau - 1 + Phi(u/h)) + h*phi(u/h)
##   psi_{tau,h}(u) = d/du rho_{tau,h}(u) = tau - 1 + Phi(u/h)          (smoothed score)
##   k_{h}(u)       = d/du psi = phi(u/h)/h                              (Hessian weight)
## As h -> 0 these converge to the check loss, its subgradient tau - 1{u<0}, and a
## point mass, recovering exact weighted quantile regression.

#' Smoothed check loss value.
#' @keywords internal
rho_smooth <- function(u, tau, h) {
  v <- u / h
  u * (tau - 1 + .pnorm(v)) + h * .dnorm(v)
}

#' Smoothed check-loss score (gradient w.r.t. the residual).
#' @keywords internal
psi_smooth <- function(u, tau, h) tau - 1 + .pnorm(u / h)

#' Smoothed-density weight (kernel at the residual), used in the Newton Hessian and
#' the sandwich variance. Read off the *fitted smoothed* density, never an ALD
#' stand-in (the f(0) lesson from mixqr).
#' @keywords internal
k_smooth <- function(u, h) .dnorm(u / h) / h

#' Dimensionless conquer-style bandwidth rate ~ scale * (log n / n)^(1/5).
#' The absolute bandwidth used in the M-step is this rate times the residual scale
#' (so the smoothing is scale-equivariant), with an absolute floor applied there.
#' @keywords internal
smooth_bandwidth <- function(n, scale = 1) {
  scale * (log(max(n, 3)) / max(n, 3))^(1 / 5)
}

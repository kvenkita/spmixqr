#' Plot a spatial quantile mixture fit
#'
#' Base-graphics maps and diagnostics. Richer \pkg{ggplot2} maps are shown in the
#' package primer.
#'
#' @param x A fitted [spmixqr] object.
#' @param which One of `"gate"` (mixing-probability map), `"coef"` (slope-surface
#'   map), `"phi"` (CAR spatial-error surface), `"class"` (classified map), or
#'   `"density"` (component error densities).
#' @param regime Which regime to map for `"gate"`/`"coef"`/`"phi"` (default 2 for
#'   gate, 1 for coef/phi).
#' @param credible For `which = "phi"`, show only the **credible** hot/cold spots
#'   (units whose `level` confidence interval for the CAR effect excludes zero); other
#'   units are greyed. Requires a fit carrying a covariance (`variance = "boot"` or
#'   `"sandwich"`).
#' @param level Confidence level used when `credible = TRUE`.
#' @param ... Passed to the underlying plotting call.
#' @return Invisibly `NULL`.
#' @export
plot.spmixqr <- function(x, which = c("gate", "coef", "phi", "class", "density"),
                         regime = NULL, credible = FALSE, level = 0.95, ...) {
  which <- match.arg(which)
  point <- !is.null(x$coords) && isTRUE(x$coords$mode == "point")
  if (which %in% c("gate", "coef", "phi", "class") && !point) {
    message("Mapping is shown for point data; for areal data use the primer's ggplot maps.")
  }
  cc <- if (point) x$coords$coords else NULL
  pal <- function(v) {
    rg <- range(v, finite = TRUE); s <- (v - rg[1]) / (diff(rg) + 1e-12)
    grDevices::rgb(s, 0.2, 1 - s)
  }
  if (which == "gate") {
    k <- if (is.null(regime)) min(2L, x$G) else regime
    v <- x$prior[, k]
    if (point) plot(cc[, 1], cc[, 2], col = pal(v), pch = 19,
                    xlab = "coord 1", ylab = "coord 2",
                    main = sprintf("Gate P(regime %d)", k), ...)
    else plot(v, type = "h", ylab = "prob", main = sprintf("Gate P(regime %d)", k))
  } else if (which == "coef") {
    if (!x$spatial_coef) { message("Flat slopes; nothing to map."); return(invisible(NULL)) }
    k <- if (is.null(regime)) 1L else regime
    ss <- slope_surface_matrix(x$coefficients, x$design, x$basis$B, 1L)[, k]
    if (point) plot(cc[, 1], cc[, 2], col = pal(ss), pch = 19,
                    xlab = "coord 1", ylab = "coord 2",
                    main = sprintf("Slope surface, regime %d", k), ...)
    else plot(ss, type = "h", main = sprintf("Slope surface, regime %d", k))
  } else if (which == "phi") {
    if (!isTRUE(x$spatial_error) || is.null(x$car)) {
      message("No CAR spatial-error term; nothing to map."); return(invisible(NULL))
    }
    k <- if (is.null(regime)) 1L else regime
    phi_u <- x$car$phi[, k]
    keep <- rep(TRUE, length(phi_u))
    main <- sprintf("CAR spatial-error surface, regime %d", k)
    if (isTRUE(credible)) {                    # show only reliable hot/cold spots
      ps <- phi_surface(x, ci = TRUE, level = level)
      ps <- ps[ps$regime == k, ]
      keep <- isTRUE_vec(ps$credible)
      main <- sprintf("Credible CAR hot/cold spots, regime %d", k)
    }
    v <- phi_u; v[!keep] <- NA                 # mask non-credible units (grey)
    vobs <- v[x$car$units$unit_idx]            # broadcast to observation rows
    if (point) {
      col <- ifelse(is.na(vobs), "grey85", pal(ifelse(is.na(vobs), 0, vobs)))
      plot(cc[, 1], cc[, 2], col = col, pch = 19,
           xlab = "coord 1", ylab = "coord 2", main = main, ...)
    } else plot(phi_u, type = "h", ylab = "phi", main = main,
                col = ifelse(keep, "black", "grey80"))
  } else if (which == "class") {
    cl <- apply(x$posterior, 1L, which.max)
    if (point) plot(cc[, 1], cc[, 2], col = cl + 1L, pch = 19,
                    xlab = "coord 1", ylab = "coord 2", main = "Classified regimes", ...)
    else plot(cl, main = "Classified regimes")
  } else {
    e <- x$y - x$fitted_q
    rng <- range(e)
    plot(stats::density(e[, 1]), xlim = rng, main = "Component error densities",
         xlab = "residual", col = 2, ...)
    for (k in seq_len(x$G)[-1]) graphics::lines(stats::density(e[, k]), col = k + 1L)
    graphics::legend("topright", legend = paste0("regime", seq_len(x$G)),
                     col = seq_len(x$G) + 1L, lty = 1, bty = "n")
  }
  invisible(NULL)
}

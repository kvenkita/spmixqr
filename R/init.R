## Multi-start initial responsibilities. The mixture likelihood is multimodal
## (Wu & Yao 2016), so spmixqr() runs several starts and keeps the best.

#' One initial n x G responsibility matrix.
#' @param y,X response and (constant) component design.
#' @param coords numeric coordinates (n x 2) or NULL/factor (areal).
#' @param G,tau number of regimes / quantile level.
#' @param type `"blocked"` (k-means on residuals + coords), `"resid"` (residuals
#'   only), or `"random"`.
#' @return An `n x G` responsibility matrix (rows sum to one).
#' @keywords internal
init_responsibilities <- function(y, X, coords, G, tau, type = "blocked") {
  n <- length(y)
  if (G == 1L) return(matrix(1, n, 1L))
  if (type == "random") {
    P <- matrix(stats::rgamma(n * G, 1), n, G)
    return(normalize_rows(P))
  }
  beta0 <- tryCatch(mixqr::weighted_rq(X, y, tau, w = rep(1, n)),
                    error = function(e) stats::coef(stats::lm.fit(X, y)))
  resid <- as.numeric(y - X %*% beta0)
  feat <- scale(resid)
  if (type == "blocked" && is.numeric(coords) && is.matrix(coords))
    feat <- cbind(feat, scale(coords))
  cl <- tryCatch(stats::kmeans(feat, centers = G, nstart = 3L)$cluster,
                 error = function(e) sample.int(G, n, replace = TRUE))
  hard <- matrix(0.2 / (G - 1), n, G)
  hard[cbind(seq_len(n), cl)] <- 0.8
  normalize_rows(hard)
}

#' A list of `nstart` initial responsibility matrices cycling through start types.
#' @keywords internal
init_starts <- function(y, X, coords, G, tau, nstart) {
  types <- c("blocked", "resid", "random")
  lapply(seq_len(nstart), function(i)
    init_responsibilities(y, X, coords, G, tau, types[((i - 1L) %% 3L) + 1L]))
}

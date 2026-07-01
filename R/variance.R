## Inference. The fast default is a classification-conditional sandwich (disclosed);
## the recommended option for spatially-dependent data is a bootstrap (xy-pairs or
## spatial-block) of the whole EM, with bootstrap labels aligned to the point
## estimate before aggregating.

#' Classification-conditional sandwich covariances (gate + per-regime coefficients).
#' @keywords internal
sandwich_vcov <- function(obj) {
  des <- obj$design
  ow <- obj$weights %||% rep(1, length(obj$y))
  wtype <- obj$weights_type %||% "sampling"
  gfit <- list(pi = obj$prior, hessian = des_gate_hessian(obj))
  Vg <- gate_sandwich_vcov(des$Z, obj$posterior, gfit, ow = ow, wtype = wtype)
  Vc <- lapply(seq_len(obj$G), function(k)
    coef_sandwich_vcov(des$Xt, obj$y, obj$tau, obj$posterior[, k],
                       obj$coefficients[, k], des$Pen_beta, obj$h,
                       ow = ow, wtype = wtype))
  list(gate = Vg, coef = Vc, method = "sandwich",
       note = "classification-conditional; use variance='boot' for reporting")
}

#' Recover the gate Hessian for the sandwich (refit at the stored gate).
#' @keywords internal
des_gate_hessian <- function(obj) {
  pen_irls_multinom(obj$design$Z, obj$posterior, obj$design$Pen_gamma,
                    w = obj$weights %||% rep(1, length(obj$y)))$hessian
}

#' Bootstrap covariance (xy-pairs or spatial-block) of the EM pipeline.
#'
#' Reuses the fixed training basis rows (resampling observations, not rebuilding the
#' basis), refits, aligns each replicate's regimes to the point estimate by matching
#' constant coefficients, and returns the empirical covariance of the stacked gate
#' and coefficient parameters. Under-covers with few blocks (Lahiri 2003); a
#' minimum-sites guard applies.
#' @keywords internal
bootstrap_vcov <- function(obj, data = NULL, B = NULL, block = NULL) {
  des <- obj$design; ctl <- obj$control %||% spmixqr_control()
  B <- if (is.null(B)) ctl$boot_B else B
  block <- if (is.null(block)) ctl$boot_block else block
  n <- length(obj$y); G <- obj$G
  kdectrl <- mixqr::mixqr_control(bandwidth = ctl$bandwidth, kde_grid = ctl$kde_grid)

  blocks <- NULL; boot_kind <- "xy-pairs"
  if (isTRUE(obj$spatial_error) && !is.null(obj$car)) {
    ## areal block bootstrap: resample contiguous blocks (connected-component /
    ## graph partition of W) rather than i.i.d. rows (which would destroy the
    ## spatial dependence the CAR term models -> anti-conservative SEs).
    blocks <- areal_blocks(obj, ctl$boot_areal_blocks)
    boot_kind <- "areal-block"
  } else if (!is.null(obj$coords) && isTRUE(obj$coords$mode == "point") && block > 1L) {
    cc <- obj$coords$coords
    nsite <- nrow(unique(cc))
    if (nsite >= ctl$min_sites) {
      bx <- cut(cc[, 1L], stats::quantile(cc[, 1L], seq(0, 1, length.out = block + 1L)),
                include.lowest = TRUE, labels = FALSE)
      by <- cut(cc[, 2L], stats::quantile(cc[, 2L], seq(0, 1, length.out = block + 1L)),
                include.lowest = TRUE, labels = FALSE)
      blocks <- interaction(bx, by, drop = TRUE)
      boot_kind <- "spatial-block"
    }
  }

  ref <- obj$beta_const
  draws <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- if (is.null(blocks)) sample.int(n, n, replace = TRUE)
           else resample_blocks(blocks)
    fit <- tryCatch(
      spatial_em_fit(obj$y[idx], des$Xt[idx, , drop = FALSE], des$Z[idx, , drop = FALSE],
                     G, obj$tau, obj$method, obj$posterior[idx, , drop = FALSE],
                     des$Pen_beta, des$Pen_gamma, obj$h, obj$spatial_coef, kdectrl, ctl,
                     spatial_error = isTRUE(obj$spatial_error)),
      error = function(e) NULL)
    if (is.null(fit)) next
    bc <- extract_const(fit$beta, des, nrow(ref))
    ord <- align_labels(bc, ref)
    fit <- permute_fit(fit, ord)
    draws[[b]] <- c(as.numeric(fit$beta), as.numeric(fit$gamma))
  }
  draws <- do.call(rbind, Filter(Negate(is.null), draws))
  if (is.null(draws) || nrow(draws) < 2L)
    return(list(gate = NULL, coef = NULL, method = "boot",
                note = "bootstrap failed to produce enough replicates"))
  Vall <- stats::cov(draws)
  ## reshape into the same layout the S3 methods expect (per-regime coef blocks +
  ## a stacked gate block), so summary()/confint() read bootstrap SEs correctly.
  P <- nrow(obj$coefficients); q1 <- nrow(obj$gamma); K <- ncol(obj$gamma)
  Vcoef <- lapply(seq_len(G), function(k) {
    idx <- ((k - 1L) * P + 1L):(k * P)
    Vall[idx, idx, drop = FALSE]
  })
  Vgate <- NULL
  if (K > 0L) {
    gidx <- (P * G + 1L):(P * G + q1 * K)
    Vgate <- Vall[gidx, gidx, drop = FALSE]
  }
  list(coef = Vcoef, gate = Vgate, all = Vall, method = "boot",
       note = sprintf("%s bootstrap, %d replicates", boot_kind, nrow(draws)))
}

#' Contiguous areal blocks for the spatial-block bootstrap.
#'
#' Partitions the CAR units into `nblk` contiguous blocks. Disconnected components are
#' kept as separate blocks; large components are further split by a coarse spatial
#' k-means on a 2-D MDS/grid embedding of the adjacency (no \pkg{igraph} dependency).
#' Returns a length-n factor (one block label per observation) so whole contiguous
#' clusters of units are resampled together (Lahiri 2003: under-covers with few
#' blocks; disclosed).
#' @keywords internal
areal_blocks <- function(obj, nblk = 8L) {
  car <- obj$car
  unit_idx <- car$units$unit_idx
  membership <- car$units$membership
  L <- nrow(car$W$W)
  comps <- sort(unique(membership))
  ublock <- integer(L)
  next_blk <- 0L
  ## embed units for spatial clustering: use degree-normalised adjacency rows
  Wm <- as.matrix(car$W$W)
  emb <- tryCatch({
    g <- Wm; d <- rowSums(g); d[d == 0] <- 1
    Lap <- diag(L) - (g / d)
    ev <- eigen(Lap, symmetric = FALSE)
    Re(ev$vectors[, order(Re(ev$values))[seq_len(min(3L, L))], drop = FALSE])
  }, error = function(e) matrix(stats::runif(L * 2L), L, 2L))
  for (cc in comps) {
    idx <- which(membership == cc)
    m <- length(idx)
    share <- max(1L, round(nblk * m / L))
    if (m <= share || share <= 1L) {
      ublock[idx] <- next_blk + 1L; next_blk <- next_blk + 1L
    } else {
      km <- tryCatch(stats::kmeans(emb[idx, , drop = FALSE], centers = share,
                                   nstart = 2L)$cluster,
                     error = function(e) rep(1L, m))
      ublock[idx] <- next_blk + km; next_blk <- next_blk + max(km)
    }
  }
  factor(ublock[unit_idx])
}

#' Resample whole spatial blocks (with replacement) to length ~ n.
#' @keywords internal
resample_blocks <- function(blocks) {
  lev <- levels(blocks)
  chosen <- sample(lev, length(lev), replace = TRUE)
  unlist(lapply(chosen, function(l) which(blocks == l)))
}

#' Align a replicate's regimes to the reference by matching constant coefficients.
#' @keywords internal
align_labels <- function(bc, ref) {
  G <- ncol(ref)
  if (G < 2L) return(1L)
  D <- matrix(0, G, G)
  for (a in seq_len(G)) for (b in seq_len(G))
    D[a, b] <- sum((bc[, a] - ref[, b])^2)
  ## greedy assignment (small G); minimise total distance
  ord <- integer(G); used <- logical(G)
  for (b in order(apply(D, 2L, min))) {
    a <- which.min(ifelse(used, Inf, D[, b]))
    ord[b] <- a; used[a] <- TRUE
  }
  ## ord[b] = which replicate-regime maps to reference b; return permutation of cols
  inv <- integer(G); inv[ord] <- seq_len(G)
  ord
}

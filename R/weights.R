## Spatial weights layer for the CAR-error module. Builds a symmetric, nonnegative
## sparse weights matrix W from sf/sp polygons (queen/rook contiguity), point
## coordinates (distance band / k-NN), an existing spdep nb/listw, or a raw matrix
## (supplied). The verified spdep idiom (Reviewer B, 2026-06-21) is reused: knn is
## DIRECTED so make.sym.nb is mandatory; style="B" (symmetric binary) is the default
## for a valid GMRF precision base (style="W" is row-standardised and asymmetric -
## lag-semantics, warned); zero.policy=TRUE everywhere for islands / empty sets.

#' Construct a spatial weights object for the CAR-error term
#'
#' Builds a symmetric, nonnegative, sparse spatial weights matrix `W` and its degree
#' matrix `D` from a variety of inputs, wrapping the result in an `spq_weights` object
#' for use as the `spatial_W` argument of [spmixqr()]. The default style `"B"`
#' (symmetric binary) gives a valid Gaussian Markov random field precision base; the
#' row-standardised style `"W"` is asymmetric (lag-semantics) and warned against for
#' the CAR penalty.
#'
#' @param x The spatial object. For `type = "queen"`/`"rook"`: an \pkg{sf}/`sp`
#'   polygon object, **or** an \pkg{spdep} `nb`/`listw` (used directly). For
#'   `type = "distance"`/`"knn"`: a two-column coordinate matrix/data frame. For
#'   `type = "supplied"`: a square nonnegative matrix (`Matrix` or base).
#' @param type Weights construction: `"queen"`/`"rook"` (polygon contiguity),
#'   `"distance"` (distance band `[d1, d2]`), `"knn"` (k nearest neighbours,
#'   symmetrised), `"supplied"` (validate a user matrix), or `"nngp"` (a sparse
#'   nearest-neighbour Gaussian-process Matern *precision* for point data: a scalable,
#'   `W`-free continuous-domain spatial-error alternative to CAR; Datta et al. 2016).
#' @param d1,d2 Lower / upper distance band (required for `"distance"`).
#' @param k Number of nearest neighbours (for `"knn"`).
#' @param style \pkg{spdep} weighting style; `"B"` (symmetric binary, default) is
#'   required for a valid GMRF precision. `"W"` (row-standardised) is asymmetric and
#'   triggers a warning (symmetrised before use).
#' @param ids Optional unit identifiers (row/column names of `W`).
#' @param m For `type = "nngp"`: number of nearest earlier-neighbours (default 10).
#' @param range For `type = "nngp"`: Matern range in raw coordinate units (`NULL` uses
#'   ~0.1 of the domain extent). Weakly identified from one realisation, so it is fixed
#'   or selected (like the smoothing penalty), not estimated (Zhang 2004).
#' @param nu For `type = "nngp"`: Matern smoothness, `0.5` (exponential) or `1.5`.
#' @param ... Currently unused.
#' @return An object of class `spq_weights`: a list with `W` (sparse symmetric
#'   `dgCMatrix`), `D` (sparse diagonal degree matrix), `ids`, `style`, `type`,
#'   `n_comp` (number of connected components, via [spdep::n.comp.nb()]), and
#'   `n_island` (number of zero-degree units). For `type = "nngp"` it additionally
#'   carries `Q` (the sparse Matern precision, used directly by the spatial-error
#'   M-step), `kind = "nngp"`, `m`, `range`, `nu`, `ordering`, and `coords`; its `W` is
#'   the symmetrised neighbour graph (used only for diagnostics and the spatial block
#'   bootstrap).
#' @references Leroux et al. (2000); \pkg{spdep} (Bivand et al.).
#' @examples
#' ## supplied matrix (no spdep input needed)
#' Wm <- matrix(0, 4, 4); Wm[1,2] <- Wm[2,1] <- Wm[2,3] <- Wm[3,2] <- 1
#' Wm[3,4] <- Wm[4,3] <- 1
#' w <- spq_weights(Wm, type = "supplied")
#' w
#' @export
spq_weights <- function(x, type = c("queen", "rook", "distance", "knn", "supplied", "nngp"),
                        d1 = 0, d2 = NULL, k = 5L, style = "B", ids = NULL,
                        m = 10L, range = NULL, nu = 0.5, ...) {
  type <- match.arg(type)
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required for spq_weights().", call. = FALSE)
  needs_spdep <- type %in% c("queen", "rook", "distance", "knn")
  if (needs_spdep && !requireNamespace("spdep", quietly = TRUE))
    stop("Package 'spdep' is required for type='", type, "'.", call. = FALSE)

  ## ---- NNGP: sparse Matern Vecchia precision for point data (early return) ----
  if (type == "nngp") {
    cc <- as.matrix(x)
    if (ncol(cc) != 2L) stop("type='nngp' needs a two-column coordinate matrix.", call. = FALSE)
    key <- apply(cc, 1L, function(z) paste(z, collapse = "_"))
    keep <- !duplicated(key)                       # merge duplicate coords to shared units
    ucc <- cc[keep, , drop = FALSE]                # distinct locations, first-occurrence order
    uid <- if (!is.null(ids)) as.character(ids)[keep] else key[keep]
    np <- nngp_precision(ucc, m = m, range = range, nu = nu)
    Q <- np$Q; W <- np$Wadj
    dimnames(Q) <- list(uid, uid); dimnames(W) <- list(uid, uid)
    return(structure(list(
      W = W, D = Matrix::Diagonal(x = as.numeric(Matrix::rowSums(W))),
      Q = Q, ids = uid, style = "nngp", type = "nngp", kind = "nngp",
      m = np$m, range = np$range, nu = np$nu, ordering = np$ordering,
      coords = ucc, n_comp = 1L, n_island = 0L),
      class = "spq_weights"))
  }

  nb <- NULL
  to_W <- function(nb, ids) {
    ## ALWAYS zero.policy=TRUE: islands / empty neighbour sets must not error.
    Wm <- spdep::nb2mat(nb, style = style, zero.policy = TRUE)   # base matrix
    W  <- Matrix::Matrix(Wm, sparse = TRUE)
    if (!is.null(ids)) dimnames(W) <- list(as.character(ids), as.character(ids))
    if (!Matrix::isSymmetric(W)) {
      if (style != "B")
        warning("style='", style, "' gives an asymmetric W; a symmetric binary ",
                "W (style='B') is required for a valid GMRF precision. ",
                "Row-standardised W is lag-semantics, not CAR. Symmetrising.",
                call. = FALSE)
      W <- Matrix::forceSymmetric((W + Matrix::t(W)) / 2)
    }
    methods::as(methods::as(W, "CsparseMatrix"), "generalMatrix")
  }

  W <- switch(type,
    queen = ,
    rook = {
      if (inherits(x, "nb"))         nb <- x
      else if (inherits(x, "listw")) nb <- x$neighbours
      else                           nb <- spdep::poly2nb(x, queen = (type == "queen"))
      ids <- if (!is.null(ids)) ids else attr(nb, "region.id")
      to_W(nb, ids)
    },
    distance = {
      cc <- as.matrix(x)
      if (is.null(d2)) stop("`distance` needs `d2` (upper band).", call. = FALSE)
      nb <- spdep::dnearneigh(cc, d1 = d1, d2 = d2)
      to_W(nb, ids)
    },
    knn = {
      cc <- as.matrix(x)
      nb <- spdep::knn2nb(spdep::knearneigh(cc, k = k))
      if (!spdep::is.symmetric.nb(nb)) nb <- spdep::make.sym.nb(nb)  # NOT optional
      to_W(nb, ids)
    },
    supplied = {
      W <- Matrix::Matrix(as.matrix(x), sparse = TRUE)
      if (nrow(W) != ncol(W)) stop("supplied W must be square.", call. = FALSE)
      if (length(W@x) && any(W@x < 0)) stop("supplied W must be nonnegative.", call. = FALSE)
      if (!Matrix::isSymmetric(W)) {
        warning("supplied W asymmetric; symmetrising (W + t(W))/2.", call. = FALSE)
        W <- Matrix::forceSymmetric((W + Matrix::t(W)) / 2)
      }
      if (!is.null(ids)) dimnames(W) <- list(as.character(ids), as.character(ids))
      methods::as(methods::as(W, "CsparseMatrix"), "generalMatrix")
    })

  d <- Matrix::rowSums(W)
  n_comp <- if (needs_spdep && !is.null(nb))
    tryCatch(spdep::n.comp.nb(nb)$nc, error = function(e) NA_integer_)
  else components_from_W(W)$nc
  structure(list(W = W, D = Matrix::Diagonal(x = as.numeric(d)),
                 ids = rownames(W), style = style, type = type,
                 n_comp = n_comp, n_island = sum(d == 0)),
            class = "spq_weights")
}

#' Connected-component count and membership from a symmetric weights matrix.
#'
#' Pure-`Matrix`/base graph traversal (no \pkg{spdep} dependency) so the `supplied`
#' path and the constraint-absorption code can find connected components of any `W`.
#' @param W A symmetric sparse/dense weights matrix.
#' @return A list with `nc` (number of components) and `membership` (length-L integer
#'   component label per unit).
#' @keywords internal
components_from_W <- function(W) {
  L <- nrow(W)
  membership <- integer(L)
  if (L == 0L) return(list(nc = 0L, membership = membership))
  ## adjacency as a list of neighbour indices
  Ws <- methods::as(W, "TsparseMatrix")
  adj <- vector("list", L)
  i <- Ws@i + 1L; j <- Ws@j + 1L
  keep <- i != j & Ws@x != 0
  i <- i[keep]; j <- j[keep]
  if (length(i)) {
    sp <- split(j, i)
    for (nm in names(sp)) adj[[as.integer(nm)]] <- sp[[nm]]
  }
  comp <- 0L
  for (start in seq_len(L)) {
    if (membership[start] != 0L) next
    comp <- comp + 1L
    stack <- start
    membership[start] <- comp
    while (length(stack)) {
      v <- stack[length(stack)]; stack <- stack[-length(stack)]
      nb <- adj[[v]]
      for (u in nb) if (membership[u] == 0L) {
        membership[u] <- comp; stack <- c(stack, u)
      }
    }
  }
  list(nc = comp, membership = membership)
}

#' Build a CAR / ICAR precision matrix from a weights object
#'
#' Forms the proper Leroux precision `Q(alpha) = alpha (D - W) + (1 - alpha) I`
#' (PSD; PD for `alpha < 1`), or the intrinsic CAR precision `D - W` plus a small
#' epsilon-ridge guard (`car = "icar"`). The proper default needs no ridge; the ICAR
#' precision is rank-deficient (null dimension = number of connected components) and
#' the ridge only makes it numerically invertible - identification of the level comes
#' from the per-component sum-to-zero constraint absorbed in [build_designs()].
#'
#' @param spqw An `spq_weights` object (or a list with `W`, `D`).
#' @param alpha Proper-CAR spatial-dependence strength in `[0, 1]` (default `0.95`).
#' @param car `"proper"` (Leroux) or `"icar"` (intrinsic).
#' @param eps Epsilon-ridge added to the ICAR precision for conditioning.
#' @return A symmetric sparse precision matrix (`dsCMatrix`).
#' @references Leroux et al. (2000).
#' @keywords internal
make_car_precision <- function(spqw, alpha = 0.95, car = c("proper", "icar"),
                               eps = 1e-6) {
  car <- match.arg(car)
  W <- spqw$W; D <- spqw$D
  Lap <- D - W
  if (car == "icar")
    return(Matrix::forceSymmetric(Lap + eps * Matrix::Diagonal(nrow(W))))
  Q <- alpha * Lap + (1 - alpha) * Matrix::Diagonal(nrow(W))      # proper: PD
  Matrix::forceSymmetric(Q)
}

#' Constraint-absorbed incidence and precision for the CAR block.
#'
#' Imposes a per-connected-component sum-to-zero constraint on `phi` (so `beta`'s
#' intercept carries the level and `phi` is the mean-zero spatial deviation), absorbing
#' the constraint into a reduced incidence basis `Rt` (`n x L'`, `L' = L - n_comp`) and
#' reduced precision `Qt` (`L' x L'`), analogous to mgcv's `absorb.cons`. A tiny
#' epsilon-ridge is added to `Qt` for conditioning. The CAR effect on the original `L`
#' units is recovered as `phi = Tmat %*% phi_reduced`, where `Tmat` (`L x L'`) is the
#' constraint null-space basis stored in the return value.
#'
#' @param R Full unit-incidence matrix `n x L` (`R[i, s(i)] = 1`), sparse.
#' @param Q Full `L x L` CAR precision from [make_car_precision()].
#' @param membership Length-L connected-component labels (from [components_from_W()]).
#' @param eps Ridge added to the reduced precision.
#' @param constrain Apply the per-component sum-to-zero (default `TRUE`, for the
#'   rank-deficient ICAR/CAR precision). `FALSE` for a proper full-rank precision (NNGP):
#'   `phi` is fit unconstrained and `Tmat` is the identity.
#' @return A list with `Rt` (constraint-absorbed incidence, sparse `n x L'`), `Qt`
#'   (reduced precision, sparse `L' x L'`), `Tmat` (`L x L'` recovery basis), and
#'   `Lp` (`L'`).
#' @keywords internal
absorb_car_constraint <- function(R, Q, membership, eps = 1e-6, constrain = TRUE) {
  L <- ncol(R)
  ## Proper full-rank precision (e.g. NNGP): no sum-to-zero. Fit phi unconstrained;
  ## the proper precision + intercept identify it (Datta et al. 2016). Tmat = identity.
  if (!isTRUE(constrain)) {
    Tmat <- methods::as(Matrix::Diagonal(L), "CsparseMatrix")
    Qt <- Matrix::forceSymmetric(Q) + eps * Matrix::Diagonal(L)
    return(list(Rt = methods::as(R, "CsparseMatrix"),
                Qt = methods::as(Qt, "CsparseMatrix"),
                Tmat = Tmat, Lp = L))
  }
  comps <- sort(unique(membership))
  ## Build the L x L' null-space basis Tmat of the per-component sum-to-zero
  ## constraint. Within each component of size m, the sum-to-zero contrast is the
  ## (m x (m-1)) basis of the orthogonal complement of the constant vector. We use
  ## the QR-based "sum-to-zero" contrast (as contr.sum-style, orthonormalised) so the
  ## reduced precision stays well-conditioned.
  cols <- vector("list", length(comps))
  rows_list <- vector("list", length(comps))
  offset <- 0L
  Tlist <- list()
  for (ci in seq_along(comps)) {
    idx <- which(membership == comps[ci])
    m <- length(idx)
    if (m <= 1L) {
      ## singleton component: a lone unit must satisfy phi = 0 (its own sum-to-zero),
      ## so it contributes no reduced column.
      next
    }
    ## orthonormal basis of {v in R^m : sum(v) = 0}: m x (m-1)
    Cmat <- diag(m) - matrix(1 / m, m, m)
    sv <- svd(Cmat)
    keep <- sv$d > 1e-8
    Bc <- sv$u[, keep, drop = FALSE]           # m x (m-1), orthonormal, columns sum to 0
    Tlist[[length(Tlist) + 1L]] <- list(idx = idx, B = Bc)
  }
  Lp <- sum(vapply(Tlist, function(t) ncol(t$B), integer(1)))
  ## assemble sparse Tmat (L x Lp)
  iT <- integer(0); jT <- integer(0); xT <- numeric(0)
  col0 <- 0L
  for (t in Tlist) {
    m <- length(t$idx); mc <- ncol(t$B)
    iT <- c(iT, rep(t$idx, times = mc))
    jT <- c(jT, rep(col0 + seq_len(mc), each = m))
    xT <- c(xT, as.numeric(t$B))
    col0 <- col0 + mc
  }
  Tmat <- if (Lp > 0L)
    Matrix::sparseMatrix(i = iT, j = jT, x = xT, dims = c(L, Lp))
  else Matrix::Matrix(0, L, 0L, sparse = TRUE)
  Rt <- R %*% Tmat
  Qt <- Matrix::forceSymmetric(Matrix::t(Tmat) %*% Q %*% Tmat)
  if (Lp > 0L) Qt <- Qt + eps * Matrix::Diagonal(Lp)
  list(Rt = methods::as(Rt, "CsparseMatrix"),
       Qt = methods::as(Qt, "CsparseMatrix"),
       Tmat = methods::as(Tmat, "CsparseMatrix"), Lp = Lp)
}

#' Unit-incidence matrix R (n x L), sparse, from an observation->unit index.
#' @param unit_idx Length-n integer index of each observation's unit (1..L).
#' @param L Number of units.
#' @return A sparse `n x L` incidence matrix.
#' @keywords internal
incidence_matrix <- function(unit_idx, L) {
  n <- length(unit_idx)
  Matrix::sparseMatrix(i = seq_len(n), j = unit_idx, x = rep(1, n), dims = c(n, L))
}

#' @export
print.spq_weights <- function(x, ...) {
  cat("<spq_weights>\n")
  if (identical(x$type, "nngp")) {
    cat(sprintf("  type = nngp (NNGP Matern precision)   nu = %s   range = %.4g\n",
                as.character(x$nu), x$range))
    cat(sprintf("  units L = %d   neighbours m = %d   precision nnz = %d\n",
                nrow(x$Q), x$m, length(x$Q@x)))
    return(invisible(x))
  }
  cat(sprintf("  type = %s   style = %s\n", x$type, x$style))
  cat(sprintf("  units L = %d   nonzero links = %d   connected components = %s\n",
              nrow(x$W), length(x$W@x), as.character(x$n_comp)))
  if (x$n_island > 0L)
    cat(sprintf("  islands (zero-degree units): %d\n", x$n_island))
  invisible(x)
}

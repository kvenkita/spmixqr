## Spatial basis + roughness penalty via mgcv smoothCon / PredictMat.
## The exact idiom verified end-to-end (build X and S; predict at new locations) for
## thin-plate, low-rank GP, and MRF smooths. Coordinates are standardised before the
## build and the transform is stored so predictions at new locations match.

#' Build a spatial basis and its roughness penalty
#'
#' Constructs a low-rank spatial basis `B(s)` (sum-to-zero constrained) and an
#' associated PSD roughness penalty `Omega` using \pkg{mgcv}. For point data a
#' thin-plate (`"tp"`) or low-rank Gaussian-process (`"gp"`) smooth; for areal data a
#' Markov-random-field (`"mrf"`) smooth from a neighbour graph.
#'
#' @param coords For point data, a two-column numeric matrix/data frame of
#'   coordinates (length n). For areal data, a length-n factor/vector of region
#'   labels (each observation's region).
#' @param areal For areal data, an \pkg{spdep} `nb` or `listw` neighbour object over
#'   the unique regions (converted to mgcv's named-list form internally). `NULL` for
#'   point data.
#' @param type Basis type; ignored (forced to `"mrf"`) when `areal` is supplied.
#' @param k Basis dimension (clamped to unique locations/regions minus one).
#' @param scale_coords Standardise point coordinates before the build.
#' @return An object of class `spmixqr_basis`: the fitted smooth `sm`, the basis
#'   matrix `B`, the penalty `Omega`, the coordinate transform, and metadata.
#' @export
spmixqr_basis <- function(coords, areal = NULL, type = c("tp", "gp", "mrf"),
                          k = 20L, scale_coords = TRUE) {
  type <- match.arg(type)
  if (!is.null(areal)) type <- "mrf"

  if (type == "mrf") {
    region <- as.factor(coords)
    nblist <- nb_to_named_list(areal, levels(region))
    nlev <- nlevels(region)
    k_use <- min(as.integer(k), nlev)
    if (k_use < as.integer(k))
      message(sprintf("k clamped to %d (number of regions).", k_use))
    df <- data.frame(.region = region)
    sm <- mgcv::smoothCon(
      mgcv::s(.region, bs = "mrf", k = k_use, xt = list(nb = nblist)),
      data = df, absorb.cons = TRUE)[[1]]
    transform <- list(type = "mrf", levels = levels(region))
  } else {
    cc <- as.matrix(coords)
    if (ncol(cc) != 2L) stop("`coords` must have two columns for point data.")
    storage.mode(cc) <- "double"
    ctr <- if (scale_coords) colMeans(cc) else c(0, 0)
    scl <- if (scale_coords) apply(cc, 2L, stats::sd) else c(1, 1)
    scl[scl == 0] <- 1
    s1 <- (cc[, 1L] - ctr[1L]) / scl[1L]
    s2 <- (cc[, 2L] - ctr[2L]) / scl[2L]
    nuniq <- nrow(unique(cbind(s1, s2)))
    k_use <- min(as.integer(k), nuniq - 1L)
    if (k_use < as.integer(k))
      message(sprintf("k clamped to %d (unique coordinate combinations - 1).", k_use))
    df <- data.frame(.s1 = s1, .s2 = s2)
    sm <- mgcv::smoothCon(
      mgcv::s(.s1, .s2, bs = type, k = k_use),
      data = df, absorb.cons = TRUE)[[1]]
    transform <- list(type = type, center = ctr, scale = scl)
  }

  Omega <- sm$S[[1]]
  Omega <- symmetrise(Omega)
  structure(list(sm = sm, B = sm$X, Omega = Omega, transform = transform,
                 type = type, r = ncol(sm$X), k = k_use),
            class = "spmixqr_basis")
}

#' Evaluate a spatial basis at new locations.
#' @param basis A `spmixqr_basis`.
#' @param newcoords New coordinates (point: two-column matrix) or region labels (mrf).
#' @return An `m x r` basis matrix aligned to the training basis columns.
#' @keywords internal
predict_basis <- function(basis, newcoords) {
  tr <- basis$transform
  if (basis$type == "mrf") {
    region <- factor(as.character(newcoords), levels = tr$levels)
    if (anyNA(region)) stop("newcoords contains regions not seen in training.")
    nd <- data.frame(.region = region)
  } else {
    cc <- as.matrix(newcoords)
    storage.mode(cc) <- "double"
    nd <- data.frame(.s1 = (cc[, 1L] - tr$center[1L]) / tr$scale[1L],
                     .s2 = (cc[, 2L] - tr$center[2L]) / tr$scale[2L])
  }
  mgcv::PredictMat(basis$sm, nd)
}

#' Convert an spdep nb/listw to mgcv's named-list neighbour form.
#' @keywords internal
nb_to_named_list <- function(areal, lev) {
  nb <- areal
  if (inherits(nb, "listw")) nb <- nb$neighbours
  if (!is.list(nb)) stop("`areal` must be an spdep nb or listw object.")
  rid <- attr(nb, "region.id")
  if (is.null(rid)) rid <- lev
  rid <- as.character(rid)
  out <- lapply(seq_along(nb), function(i) {
    nn <- nb[[i]]
    nn <- nn[nn > 0]                       # spdep codes "no neighbours" as 0
    as.character(rid[nn])
  })
  names(out) <- rid
  if (!all(lev %in% names(out)))
    stop("Neighbour object's region ids do not match the region labels in `coords`. ",
         "Ensure the nb/listw `region.id` equals levels(factor(coords)).",
         call. = FALSE)
  out[lev]                                  # order to match the factor levels
}

#' @export
print.spmixqr_basis <- function(x, ...) {
  cat(sprintf("<spmixqr_basis> type=%s  dimension r=%d (k=%d)\n",
              x$type, x$r, x$k))
  invisible(x)
}

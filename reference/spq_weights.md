# Construct a spatial weights object for the CAR-error term

Builds a symmetric, nonnegative, sparse spatial weights matrix `W` and
its degree matrix `D` from a variety of inputs, wrapping the result in
an `spq_weights` object for use as the `spatial_W` argument of
[`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md).
The default style `"B"` (symmetric binary) gives a valid Gaussian Markov
random field precision base; the row-standardised style `"W"` is
asymmetric (lag-semantics) and warned against for the CAR penalty.

## Usage

``` r
spq_weights(
  x,
  type = c("queen", "rook", "distance", "knn", "supplied", "nngp"),
  d1 = 0,
  d2 = NULL,
  k = 5L,
  style = "B",
  ids = NULL,
  m = 10L,
  range = NULL,
  nu = 0.5,
  ...
)
```

## Arguments

- x:

  The spatial object. For `type = "queen"`/`"rook"`: an sf/`sp` polygon
  object, **or** an spdep `nb`/`listw` (used directly). For
  `type = "distance"`/`"knn"`: a two-column coordinate matrix/data
  frame. For `type = "supplied"`: a square nonnegative matrix (`Matrix`
  or base).

- type:

  Weights construction: `"queen"`/`"rook"` (polygon contiguity),
  `"distance"` (distance band `[d1, d2]`), `"knn"` (k nearest
  neighbours, symmetrised), `"supplied"` (validate a user matrix), or
  `"nngp"` (a sparse nearest-neighbour Gaussian-process Matern
  *precision* for point data: a scalable, `W`-free continuous-domain
  spatial-error alternative to CAR; Datta et al. 2016).

- d1, d2:

  Lower / upper distance band (required for `"distance"`).

- k:

  Number of nearest neighbours (for `"knn"`).

- style:

  spdep weighting style; `"B"` (symmetric binary, default) is required
  for a valid GMRF precision. `"W"` (row-standardised) is asymmetric and
  triggers a warning (symmetrised before use).

- ids:

  Optional unit identifiers (row/column names of `W`).

- m:

  For `type = "nngp"`: number of nearest earlier-neighbours (default
  10).

- range:

  For `type = "nngp"`: Matern range in raw coordinate units (`NULL` uses
  ~0.1 of the domain extent). Weakly identified from one realisation, so
  it is fixed or selected (like the smoothing penalty), not estimated
  (Zhang 2004).

- nu:

  For `type = "nngp"`: Matern smoothness, `0.5` (exponential) or `1.5`.

- ...:

  Currently unused.

## Value

An object of class `spq_weights`: a list with `W` (sparse symmetric
`dgCMatrix`), `D` (sparse diagonal degree matrix), `ids`, `style`,
`type`, `n_comp` (number of connected components, via
[`spdep::n.comp.nb()`](https://r-spatial.github.io/spdep/reference/compon.html)),
and `n_island` (number of zero-degree units). For `type = "nngp"` it
additionally carries `Q` (the sparse Matern precision, used directly by
the spatial-error M-step), `kind = "nngp"`, `m`, `range`, `nu`,
`ordering`, and `coords`; its `W` is the symmetrised neighbour graph
(used only for diagnostics and the spatial block bootstrap).

## References

Leroux et al. (2000); spdep (Bivand et al.).

## Examples

``` r
## supplied matrix (no spdep input needed)
Wm <- matrix(0, 4, 4); Wm[1,2] <- Wm[2,1] <- Wm[2,3] <- Wm[3,2] <- 1
Wm[3,4] <- Wm[4,3] <- 1
w <- spq_weights(Wm, type = "supplied")
w
#> <spq_weights>
#>   type = supplied   style = B
#>   units L = 4   nonzero links = 6   connected components = 1
```

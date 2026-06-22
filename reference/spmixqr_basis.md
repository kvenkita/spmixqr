# Build a spatial basis and its roughness penalty

Constructs a low-rank spatial basis `B(s)` (sum-to-zero constrained) and
an associated PSD roughness penalty `Omega` using mgcv. For point data a
thin-plate (`"tp"`) or low-rank Gaussian-process (`"gp"`) smooth; for
areal data a Markov-random-field (`"mrf"`) smooth from a neighbour
graph.

## Usage

``` r
spmixqr_basis(
  coords,
  areal = NULL,
  type = c("tp", "gp", "mrf"),
  k = 20L,
  scale_coords = TRUE
)
```

## Arguments

- coords:

  For point data, a two-column numeric matrix/data frame of coordinates
  (length n). For areal data, a length-n factor/vector of region labels
  (each observation's region).

- areal:

  For areal data, an spdep `nb` or `listw` neighbour object over the
  unique regions (converted to mgcv's named-list form internally).
  `NULL` for point data.

- type:

  Basis type; ignored (forced to `"mrf"`) when `areal` is supplied.

- k:

  Basis dimension (clamped to unique locations/regions minus one).

- scale_coords:

  Standardise point coordinates before the build.

## Value

An object of class `spmixqr_basis`: the fitted smooth `sm`, the basis
matrix `B`, the penalty `Omega`, the coordinate transform, and metadata.

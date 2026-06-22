# NNGP sparse Matern precision and its neighbour graph.

Builds `Q = (I - B)' D^{-1} (I - B)` (Datta et al. 2016) on the supplied
point coordinates, returning the sparse precision (in the input row
order) plus the symmetrised neighbour adjacency (used by diagnostics /
the spatial-block bootstrap).

## Usage

``` r
nngp_precision(coords, m = 10L, range = NULL, nu = 0.5)
```

## Arguments

- coords:

  An `n x 2` coordinate matrix (distinct locations; raw units).

- m:

  Number of nearest earlier-neighbours (default 10).

- range:

  Matern range in raw coordinate units; `NULL` =\> ~0.1 \* domain
  extent.

- nu:

  Matern smoothness in `{0.5, 1.5}`.

## Value

A list with `Q` (sparse symmetric `dgCMatrix`, input order), `Wadj`
(sparse symmetric 0/1 neighbour adjacency), `range`, `m`, `nu`,
`ordering` (max-min order).

## References

Datta, Banerjee, Finley & Gelfand (2016, JASA); Guinness (2018,
Technometrics); Zhang (2004, JASA).

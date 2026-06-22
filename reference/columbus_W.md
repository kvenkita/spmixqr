# Queen-contiguity weights for the Columbus neighbourhoods

A sparse, symmetric binary (`style = "B"`) queen-contiguity weights
matrix over the 49 Columbus neighbourhoods of
[columbus_crime](https://kvenkita.github.io/spmixqr/reference/columbus_crime.md)
(118 links, one connected component). Pass to
[`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
via `spatial_W = spq_weights(columbus_W, type = "supplied")`.

## Usage

``` r
columbus_W
```

## Format

A 49 x 49 sparse `dgCMatrix` (118 nonzero links).

## Source

Queen contiguity (`spdep::poly2nb(queen = TRUE)`) of the spData Columbus
shapefile; built offline in `data-raw/make_columbus.R`.

## Examples

``` r
data(columbus_W)
dim(columbus_W)
#> Loading required package: Matrix
#> NULL
```

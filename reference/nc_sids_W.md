# Queen-contiguity weights for the North Carolina SIDS counties

A sparse, symmetric binary (`style = "B"`) queen-contiguity weights
matrix over the 100 North Carolina counties of
[nc_sids](https://kvenkita.github.io/spmixqr/reference/nc_sids.md), with
row/column names equal to `nc_sids$county`. Pass to
[`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
via `spatial_W = spq_weights(nc_sids_W, type = "supplied")` (or directly
as `spatial_W`).

## Usage

``` r
nc_sids_W
```

## Format

A 100 x 100 sparse `dgCMatrix` (490 nonzero links, one connected
component, no islands).

## Source

Queen contiguity (`spdep::poly2nb(queen = TRUE)`) of the sf North
Carolina shapefile; built offline in `data-raw/make_nc_tracts.R`.

## Examples

``` r
data(nc_sids_W)
dim(nc_sids_W)
#> [1] 100 100
```

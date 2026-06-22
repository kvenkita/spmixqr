# North Carolina SIDS counts (areal census example)

Sudden-infant-death (SIDS) counts and live births for the 100 counties
of North Carolina, 1974–78, the canonical areal/lattice spatial dataset.
Shipped as a tidy data frame with county centroids and a derived log
SIDS rate, to illustrate the CAR spatial-error quantile regression on
contiguous (polygon) areal units. The matching queen-contiguity weights
are in
[nc_sids_W](https://kvenkita.github.io/spmixqr/reference/nc_sids_W.md).

## Usage

``` r
nc_sids
```

## Format

A data frame with 100 rows (counties) and 9 variables:

- county:

  County name (the spatial unit id; matches `nc_sids_W` row/col names).

- births:

  Live births, 1974–78 (`BIR74`).

- sid:

  SIDS deaths, 1974–78 (`SID74`).

- nwbirths:

  Non-white live births, 1974–78 (`NWBIR74`).

- east, north:

  County-centroid longitude / latitude.

- sids_rate:

  SIDS per 1000 births (`1000 * (sid + 0.5) / (births + 1)`).

- log_sids:

  `log(sids_rate)` – the modelled response.

- pnw:

  Proportion of non-white births (`nwbirths / (births + 1)`).

## Source

The `nc.shp` North Carolina SIDS shapefile distributed with the sf
package (originally Cressie & Read 1985; Cressie 1993, *Statistics for
Spatial Data*). Built offline from
`sf::st_read(system.file("shape/nc.shp", package = "sf"))` and
redistributed as a built-in for reproducible examples; see
`data-raw/make_nc_tracts.R`.

## Examples

``` r
data(nc_sids)
data(nc_sids_W)
summary(nc_sids$log_sids)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -0.9813  0.3562  0.7394  0.7198  1.0610  2.2891 
```

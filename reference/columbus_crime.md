# Columbus, Ohio neighbourhood crime

Residential crime in 49 neighbourhoods of Columbus, Ohio (1980), the
canonical spatial-econometrics dataset (Anselin 1988). Shipped as a tidy
data frame with neighbourhood centroids, to illustrate spatial-error
(CAR) quantile regression on contiguous areal units. The matching
queen-contiguity weights are in
[columbus_W](https://kvenkita.github.io/spmixqr/reference/columbus_W.md).

## Usage

``` r
columbus_crime
```

## Format

A data frame with 49 rows (neighbourhoods) and 6 variables:

- id:

  Neighbourhood id (the spatial unit; matches `columbus_W` row/col
  order).

- crime:

  Residential burglaries and vehicle thefts per 1000 households.

- income:

  Household income (USD 1000).

- hoval:

  Housing value (USD 1000).

- x, y:

  Neighbourhood-polygon centroids (planar map units).

## Source

The `columbus` data distributed with the spData package (Anselin, L.
1988, *Spatial Econometrics: Methods and Models*). Built offline from
the spData `columbus.gpkg` shapefile and redistributed as a built-in for
reproducible examples; see `data-raw/make_columbus.R`.

## Examples

``` r
data(columbus_crime)
summary(columbus_crime$crime)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.1783 20.0485 34.0008 35.1288 48.5855 68.8920 
```

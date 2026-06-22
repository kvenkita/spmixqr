# Apply stored Spatial+ smooths to residualise new covariates (predict stage).

Mirrors
[`spatial_plus_residualize()`](https://kvenkita.github.io/spmixqr/reference/spatial_plus_residualize.md)
on `newdata`: subtracts each stored smooth's prediction (at the new
coordinates / regions) from the corresponding covariate column.

## Usage

``` r
spatial_plus_apply(Xnew, smooths, geo_new)
```

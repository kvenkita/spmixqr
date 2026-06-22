# Residualise covariates against a spatial smooth (Spatial+ stage 1).

For each non-intercept column of the component design, regress it on a
spatial smooth (a thin-plate smooth of point coordinates, or a
Markov-random-field smooth over areal regions) and replace it by the
residual. Returns the residualised design, the fitted smooths (needed to
residualise `newdata` at prediction), and the per-covariate spatial
R-squared removed.

## Usage

``` r
spatial_plus_residualize(X, slope_idx, geo, k = NULL)
```

## Arguments

- X:

  The component design matrix (with intercept).

- slope_idx:

  Integer indices of the non-intercept columns to residualise.

- geo:

  The resolved geography
  ([`resolve_coords()`](https://kvenkita.github.io/spmixqr/reference/resolve_coords.md)
  output): point `coords` or areal `region` + `areal` nb.

- k:

  Basis dimension for the smooth (`NULL` = a generous default richer
  than the spatial-error resolution).

## Value

A list with `X` (residualised design), `smooths` (named list of fitted
`gam`s, one per residualised covariate; `NULL` for skipped constant
columns), and `spatialR2` (a data frame of the spatial R-squared removed
per covariate).

## References

Dupont, Wood & Augustin (2022, Biometrics).

# Spatial coefficient surfaces

Evaluates each regime's covariate-effect slope surface(s) `beta_kj(s)`
over a set of locations (the training locations by default), returning a
tidy data frame for mapping.

## Usage

``` r
coef_surface(object, newcoords = NULL, covariate = 1L)
```

## Arguments

- object:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object.

- newcoords:

  Optional coordinates (point) or region labels (areal); defaults to the
  training locations.

- covariate:

  Which slope covariate (index among the non-intercept terms).

## Value

A data frame with coordinates, `regime`, and `slope`.

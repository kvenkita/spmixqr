# Predictions from a spatial quantile mixture model

Predictions from a spatial quantile mixture model

## Usage

``` r
# S3 method for class 'spmixqr'
predict(
  object,
  newdata = NULL,
  newcoords = NULL,
  type = c("prob", "class", "posterior", "quantile", "quantile_byclass"),
  ...
)
```

## Arguments

- object:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object.

- newdata:

  Optional data frame of predictors. If `NULL`, in-sample.

- newcoords:

  Coordinates (point: two-column) or region labels (areal) for
  `newdata`. Required when `newdata` is supplied and the model is
  spatial.

- type:

  One of `"prob"` (gate probabilities over space, `n x G`), `"class"`
  (most probable regime by the gate), `"posterior"` (responsibilities;
  needs the response in `newdata`), `"quantile"` (gate-mixed conditional
  quantile), or `"quantile_byclass"` (per-regime conditional quantile,
  `n x G`).

- ...:

  Ignored.

## Value

A vector or matrix depending on `type`.

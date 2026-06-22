# Spatial gate surfaces

Evaluates the mixing probabilities over space, returning a tidy data
frame.

## Usage

``` r
gate_surface(object, newcoords = NULL, newdata = NULL)
```

## Arguments

- object:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object.

- newcoords:

  Optional coordinates / region labels; defaults to training.

- newdata:

  Optional data frame for gating covariates (if any).

## Value

A data frame with coordinates, `regime`, and `prob`.

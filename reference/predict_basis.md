# Evaluate a spatial basis at new locations.

Evaluate a spatial basis at new locations.

## Usage

``` r
predict_basis(basis, newcoords)
```

## Arguments

- basis:

  A `spmixqr_basis`.

- newcoords:

  New coordinates (point: two-column matrix) or region labels (mrf).

## Value

An `m x r` basis matrix aligned to the training basis columns.

# Symmetrise a matrix (numerical hygiene for covariances).

Coerces a Matrix S4 object to a base dense matrix first, so the result
is always a plain numeric matrix (covariances downstream are read with
base [`diag()`](https://rdrr.io/r/base/diag.html)).

## Usage

``` r
symmetrise(A)
```

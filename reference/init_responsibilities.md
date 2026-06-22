# One initial n x G responsibility matrix.

One initial n x G responsibility matrix.

## Usage

``` r
init_responsibilities(y, X, coords, G, tau, type = "blocked")
```

## Arguments

- y, X:

  response and (constant) component design.

- coords:

  numeric coordinates (n x 2) or NULL/factor (areal).

- G, tau:

  number of regimes / quantile level.

- type:

  `"blocked"` (k-means on residuals + coords), `"resid"` (residuals
  only), or `"random"`.

## Value

An `n x G` responsibility matrix (rows sum to one).

# Global component order from a constant-coefficient matrix.

Global component order from a constant-coefficient matrix.

## Usage

``` r
order_components(beta_const, label_order = "slope", order_var = 1L)
```

## Arguments

- beta_const:

  `p x G` constant coefficients (intercept + constant slopes).

- label_order:

  `"slope"` or `"intercept"`.

- order_var:

  Slope index used when ordering by slope (1 = first slope).

## Value

An integer permutation (ascending key).

# Constant coefficient matrix (p x G) from the augmented coefficients.

Returns only the beta (intercept + constant-slope) rows; the
spatial-slope basis columns and the CAR block (`des$car_block`) are
excluded, so relabelling, ordering, and label alignment never see the
high-dimensional CAR noise.

## Usage

``` r
extract_const(beta, des, p)
```

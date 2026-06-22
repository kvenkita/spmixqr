# Effective degrees of freedom (component smoothers + gate smoother).

The component trace `tr[(X'WpX + P)^{-1} X'WpX]` is computed sparsely
when the CAR block is present (the design is a sparse `Matrix`), adding
the CAR-smoother trace.

## Usage

``` r
compute_edf(fit, Xt, Z, Pen_beta, Pen_gamma, G, des, spatial_coef)
```

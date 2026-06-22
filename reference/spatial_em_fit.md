# One penalised-EM fit from an initial responsibility matrix.

When `spatial_error = TRUE` the component M-step is **always** routed
through the penalised CAR solver
([`pen_smooth_wqr_car()`](https://kvenkita.github.io/spmixqr/reference/pen_smooth_wqr_car.md))
on the augmented `[X, Rt]` design with the `lambda_phi * Qt` penalty
block, regardless of `spatial_coef` (the unpenalised `weighted_rq`
branch would silently drop the CAR penalty for the `G = 1` /
constant-slope paths).

## Usage

``` r
spatial_em_fit(
  y,
  Xt,
  Z,
  G,
  tau,
  method,
  p_init,
  Pen_beta,
  Pen_gamma,
  h,
  spatial_coef,
  kdectrl,
  control,
  spatial_error = FALSE
)
```

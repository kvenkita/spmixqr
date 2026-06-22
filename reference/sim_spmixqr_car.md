# Simulate from a CAR spatial-error quantile mixture on a lattice.

Lattice (`L = lattice^2` units, rook contiguity); each regime carries a
proper-CAR spatial random effect drawn from `N(0, [Q(alpha)]^{-1})`,
scaled to SD `car_rho` and centred (mean-zero), plus a constant slope on
`x`. The conditional tau-quantile is exact by construction. Returns
`spatial_W`, `region`, `truth$phi`, `truth$unit`.

## Usage

``` r
sim_spmixqr_car(
  n,
  G,
  tau,
  coef_slope,
  sep,
  sd,
  error,
  lattice,
  car_alpha,
  car_rho
)
```

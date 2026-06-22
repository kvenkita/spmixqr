# Assemble augmented component / gate designs and their penalties.

When `car` is non-`NULL` (the CAR spatial-error path) it must be a list
with the constraint-absorbed incidence `Rt` (`n x L'`), reduced
precision `Qt` (`L' x L'`), and `lambda` (the CAR penalty weight
`lambda_phi`). The CAR columns are appended to the component design `Xt`
and the penalty `lambda * Qt` is placed on the CAR block; `car_block`
records the CAR column indices (so relabelling / ordering ignore them).

## Usage

``` r
build_designs(
  X,
  W,
  B,
  Omega,
  slope_idx,
  spatial_gate,
  spatial_coef,
  lam_g,
  lam_b,
  gate_ridge,
  r,
  car = NULL
)
```

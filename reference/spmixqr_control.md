# Control parameters for [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)

Tuning knobs for the penalised EM, the spatial basis, label handling,
and the bootstrap. All have sensible defaults.

## Usage

``` r
spmixqr_control(
  nstart = 5L,
  maxit = 200L,
  tol = 1e-05,
  basis_type = c("tp", "gp", "mrf"),
  k = 20L,
  scale_coords = TRUE,
  label_order = c("slope", "intercept"),
  order_var = 1L,
  gate_ridge = 0.001,
  gate_maxit = 50L,
  gate_tol = 1e-08,
  sm_scale = 1,
  sm_floor = 0.001,
  coef_maxit = 50L,
  coef_tol = 1e-07,
  bandwidth = NULL,
  kde_grid = 512L,
  boot_B = 200L,
  boot_block = 4L,
  boot_areal_blocks = 8L,
  min_sites = 12L,
  lambda_error_default = 1,
  car_alpha_grid = c(0.5, 0.9, 0.95, 0.99),
  trace = FALSE,
  seed = NULL
)
```

## Arguments

- nstart:

  Number of EM starts (the mixture likelihood is multimodal; Wu & Yao
  2016 stress initial-value sensitivity). The best start by penalised
  smoothed objective is kept.

- maxit, tol:

  EM iteration cap and relative-change tolerance.

- basis_type:

  Spatial basis: `"tp"` (thin-plate), `"gp"` (low-rank Gaussian process)
  for point data, or `"mrf"` (Markov random field) for areal data.

- k:

  Spatial-basis dimension (clamped to the number of unique
  locations/regions minus one).

- scale_coords:

  Standardise coordinates before building the basis (kept as a stored
  transform for prediction at new locations). Recommended `TRUE`.

- label_order:

  Component ordering key for relabelling: `"slope"` (default, per-regime
  average ordering-covariate slope) or `"intercept"`.

- order_var:

  Index (among the slope covariates) used as the ordering covariate;
  default `1` (the first non-intercept covariate).

- gate_ridge:

  Ridge on the non-spatial gate coefficients (also stabilises the
  ICAR/MRF null space). Matches the `mixqrgate` default for the
  no-spatial reduction.

- gate_maxit, gate_tol:

  Inner gate Newton/IRLS controls.

- sm_scale, sm_floor:

  Smoothing-bandwidth tuning constant and floor for the
  convolution-smoothed component M-step.

- coef_maxit, coef_tol:

  Inner penalised-Newton controls for the component step.

- bandwidth, kde_grid:

  Passed to
  [`mixqr::mixqr_control()`](https://kvenkita.github.io/mixqr/reference/mixqr_control.html)
  for the KDE density path.

- boot_B, boot_block:

  Bootstrap replicates and number of spatial blocks (per axis) for the
  spatial-block bootstrap.

- boot_areal_blocks:

  Target number of contiguous blocks for the areal (CAR) block bootstrap
  (connected-component / graph partition of `W`).

- min_sites:

  Minimum distinct locations required before a block bootstrap is
  attempted (it under-covers with few blocks; Lahiri 2003).

- lambda_error_default:

  Default CAR penalty `lambda_phi` when `lambda_error` is `NULL` (the
  CAR effect's roughness; select by BIC with
  [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md)).

- car_alpha_grid:

  Candidate proper-CAR `alpha` values for the optional coarse profile in
  [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md)
  (alpha is fixed by default, not check-loss selected).

- trace:

  Logical; print EM progress.

- seed:

  Optional RNG seed (honoured throughout for reproducibility).

## Value

A list of control parameters.

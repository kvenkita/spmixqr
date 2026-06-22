# Select regimes and smoothing parameters

Fits [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
over a grid of `G` and/or roughness penalties and returns the best by
BIC (post-hoc unsmoothed log-likelihood plus effective df) or by K-fold
held-out predictive check loss. For the KDE density path the
log-likelihood is unavailable, so `criterion = "cv"` (the check-loss
surrogate) is used.

## Usage

``` r
spmixqr_select(
  formula,
  data,
  coords = NULL,
  areal = NULL,
  tau = 0.5,
  gating = ~1,
  spatial_gate = TRUE,
  spatial_coef = TRUE,
  spatial_error = FALSE,
  spatial_W = NULL,
  car = c("proper", "icar"),
  car_alpha = 0.95,
  spatial_plus = FALSE,
  spatial_plus_k = NULL,
  method = c("ald", "kde"),
  G_grid = 2:3,
  lambda_gate_grid = c(0.1, 1, 10),
  lambda_coef_grid = c(0.1, 1, 10),
  lambda_error_grid = c(0.1, 1, 10),
  criterion = c("bic", "cv"),
  folds = 5L,
  control = spmixqr_control()
)
```

## Arguments

- formula, data, coords, areal, tau, gating, spatial_gate, spatial_coef,
  method:

  Passed to
  [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md).

- spatial_error, spatial_W, car, car_alpha:

  CAR spatial-error settings passed to
  [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md).
  When `spatial_error = TRUE` the `lambda_error_grid` is searched (BIC,
  with the CAR-effective-df term) and `spatial_gate` is forced off (the
  guardrail); `car_alpha` is fixed (not selected by the check loss).

- spatial_plus, spatial_plus_k:

  Spatial+ confounding-safeguard settings passed to
  [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  (the residualisation is refit within each CV training fold, so there
  is no leakage). NNGP `range` selection is not searched here (use
  `criterion = "bic"` comparing manually-built NNGP weights objects).

- G_grid:

  Candidate regime counts.

- lambda_gate_grid, lambda_coef_grid, lambda_error_grid:

  Candidate penalties.

- criterion:

  `"bic"` or `"cv"`.

- folds:

  Number of CV folds (for `criterion = "cv"`).

- control:

  A
  [`spmixqr_control()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_control.md)
  list.

## Value

A list with the best fit, the chosen settings, and the score table.

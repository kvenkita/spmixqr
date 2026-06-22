# Fit a spatial finite mixture of quantile regressions

Fits `G` latent regimes, each a `tau`-th quantile regression, whose
mixing probabilities vary over space through a spatial basis and
(optionally) whose covariate-effect slopes are spatially varying
surfaces. Estimation is a penalised EM. With `G = 1` this is a penalised
spatially-varying-coefficient quantile regression. See the package
primer for a worked example.

## Usage

``` r
spmixqr(
  formula,
  data,
  coords = NULL,
  areal = NULL,
  G = 2L,
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
  lambda_gate = NULL,
  lambda_coef = NULL,
  lambda_error = NULL,
  variance = c("sandwich", "boot", "none"),
  basis = NULL,
  control = spmixqr_control()
)
```

## Arguments

- formula:

  Component model formula, e.g. `y ~ x1 + x2`.

- data:

  A data frame.

- coords:

  Point coordinates: a two-column matrix/data frame, or the names of two
  columns in `data`. For areal data, the region label per row (a factor
  or a column name), paired with `areal`.

- areal:

  An spdep `nb`/`listw` neighbour object over the regions (areal data).
  `NULL` for point data.

- G:

  Number of regimes.

- tau:

  Quantile level in (0, 1).

- gating:

  Gating-covariate formula (besides space), e.g. `~ z`. Default `~ 1`.

- spatial_gate:

  Logical; let the mixing probabilities vary over space.

- spatial_coef:

  Logical; let component slopes vary over space (scalar intercepts; see
  Details).

- spatial_error:

  Logical; add a per-regime conditional-autoregressive (CAR) spatial
  random effect `phi_k` (a mean-zero spatial-level surface) to each
  component. When `TRUE` the spatial gate is turned off by default (a
  free spatial gate aliases with the CAR intercept surface; see
  Details); an explicit `spatial_gate = TRUE` then raises a guardrail
  error.

- spatial_W:

  A
  [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  object or a square weights matrix for the CAR term. `NULL`
  auto-builds: queen contiguity from `areal` (an `nb`/`listw`), else a
  k-nearest-neighbour graph from point `coords`.

- car:

  CAR precision family: `"proper"` (Leroux `alpha(D-W)+(1-alpha)I`,
  default) or `"icar"` (intrinsic `D-W` with a sum-to-zero constraint).

- car_alpha:

  Proper-CAR spatial-dependence strength in `[0, 1]` (default `0.95`).
  Fixed by default (weakly identified; disclosed in
  [`summary()`](https://rdrr.io/r/base/summary.html)), not selected by
  the check loss.

- spatial_plus:

  Logical; apply the **Spatial+** confounding safeguard (Dupont, Wood &
  Augustin 2022): residualise each covariate against a spatial smooth
  and fit on the residuals, so a smoothly-spatial covariate no longer
  competes with the spatial random effect. The reported slopes are then
  the effect of the *non-spatial* part of each covariate. Deconfounds
  only to the extent the residualisation smooth out-resolves the
  (penalised) spatial term (Frisch–Waugh–Lovell); the default smooth is
  made richer for this reason. Composes with `spatial_error` (CAR or
  NNGP).

- spatial_plus_k:

  Basis dimension for the Spatial+ covariate smooth (`NULL` uses a
  generous default richer than the spatial-error resolution).

- method:

  Component density: `"ald"` (asymmetric Laplace) or `"kde"` (Wu & Yao
  constrained kernel density).

- lambda_gate, lambda_coef, lambda_error:

  Roughness / CAR penalties for the gate, slope surfaces, and CAR effect
  (`lambda_error` = `lambda_phi`). `NULL` uses the control default;
  choose with
  [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md).

- variance:

  Inference: `"sandwich"` (fast, classification-conditional; default),
  `"boot"` (spatial-block/xy bootstrap; recommended for reporting), or
  `"none"`.

- basis:

  Advanced override: a prebuilt
  [`spmixqr_basis()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_basis.md).
  Usually `NULL` (built internally from `coords`/`areal`).

- control:

  A
  [`spmixqr_control()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_control.md)
  list.

## Value

An object of class `spmixqr`.

## Details

A free spatial intercept surface is not identified separately from the
spatial gate (both move the marginal quantile level across space). v1
therefore models spatial level/membership through the gate and spatial
covariate effects through the slope surfaces, with scalar per-regime
intercepts. Component identities are ordered by a sign-stable separating
slope; the `label_stability` diagnostic warns when slope surfaces cross
(a single global label is then only approximately coherent).

## References

Reich, Fuentes & Dunson (2011); Wu & Yao (2016); Fernandes, Guerre &
Horta (2021).

## Examples

``` r
set.seed(1)
d <- sim_spmixqr(n = 200, G = 2, tau = 0.5)
fit <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
               variance = "none", control = spmixqr_control(nstart = 2L))
fit
#> Spatial mixture of quantile regressions (spmixqr)
#>   G = 2 regimes,  tau = 0.5,  method = ald
#>   spatial gate: TRUE   spatial slopes: TRUE   basis: tp (r=19)
#> 
#> Constant component coefficients (intercept + average slopes):
#>             regime1 regime2
#> (Intercept)  2.9676 -0.1565
#> x           -2.5941  1.2814
#> 
#>   logLik = -353.24   edf = 22.10   AIC = 750.7   BIC = 823.6
```

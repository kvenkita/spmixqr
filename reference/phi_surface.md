# CAR spatial-error surfaces (with uncertainty)

Returns the per-regime CAR spatial random effect `phi_k` (the mean-zero
spatial-level deviation, on the response scale at the fitted quantile)
as a tidy data frame for mapping, mirroring
[`coef_surface()`](https://kvenkita.github.io/spmixqr/reference/coef_surface.md)
and
[`gate_surface()`](https://kvenkita.github.io/spmixqr/reference/gate_surface.md).
One row per (unit, regime). With `ci = TRUE` it adds standard errors and
a confidence interval, and flags the `credible` units whose interval
excludes zero (reliable hot/cold spots, as distinct from smoothing).
With `scale = "exp"` it adds the multiplicative deviation
`mult = exp(phi)` (interpretable for log-outcome models: e.g. `mult`
1.18 means the outcome runs about 18% above what the covariates predict
there).

## Usage

``` r
phi_surface(
  object,
  newunits = NULL,
  ci = FALSE,
  level = 0.95,
  scale = c("link", "exp")
)
```

## Arguments

- object:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object with `spatial_error = TRUE`.

- newunits:

  Optional unit identifiers to restrict / reorder the output (validated
  against the CAR unit ids). `NULL` uses all training units.

- ci:

  Add `se`, `lower`, `upper`, and `credible` columns. Requires the fit
  to carry a covariance (`variance = "boot"`, recommended, or
  `"sandwich"`).

- level:

  Confidence level for the interval.

- scale:

  `"link"` (the response/quantile scale, default) or `"exp"` (add the
  multiplicative deviation `exp(phi)`, for log-outcome models).

## Value

A data frame with `unit`, `regime`, `phi`, and (per `ci`/`scale`) `se`,
`lower`, `upper`, `credible`, `mult`, `mult_lower`, `mult_upper`.

## Examples

``` r
# \donttest{
d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE, lattice = 6,
                 seed = 1)
fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
               spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
               control = spmixqr_control(nstart = 1L, seed = 1))
#> spatial_error = TRUE: turning the spatial gate off (the free spatial gate aliases with the CAR intercept surface). Set gating = ~z for a covariate gate.
head(phi_surface(fit, ci = TRUE))
#>   unit regime         phi        se       lower       upper credible
#> 1    1      1 -0.42836396 0.2056924 -0.83151370 -0.02521422     TRUE
#> 2    2      1  0.19297118 0.2737445 -0.34355819  0.72950054    FALSE
#> 3    3      1  0.34424010 0.2111965 -0.06969738  0.75817759    FALSE
#> 4    4      1  1.10154957 0.2532079  0.60527121  1.59782793     TRUE
#> 5    5      1  0.82522616 0.1792360  0.47393004  1.17652228     TRUE
#> 6    6      1 -0.01626368 0.5369815 -1.06872817  1.03620081    FALSE
# }
```

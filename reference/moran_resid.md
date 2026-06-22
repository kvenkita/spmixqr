# Permutation Moran's I on the spatial-unit residuals of an spmixqr fit

Computes Moran's I of the responsibility-weighted quantile residual
aggregated to the spatial unit, with a permutation p-value (random
relabelling of units). For a CAR fit (`spatial_error = TRUE`) it reports
the residual Moran's I *after* the CAR term and, by re-fitting without
it, the value *before* - the negative-control / power diagnostic. For a
non-spatial-error fit it reports the single residual Moran's I against
the supplied weights.

## Usage

``` r
moran_resid(object, spatial_W = NULL, nsim = 999L)
```

## Arguments

- object:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object.

- spatial_W:

  A
  [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  object or weights matrix. Defaults to the fit's CAR weights
  (`object$car$W`) when present; otherwise required.

- nsim:

  Number of permutations for the p-value.

## Value

A list of class `spq_moran` with `before`/`after` (each a list with
`statistic`, `p_value`, `n_units`), or a single `statistic`/`p_value`
for a non-CAR fit.

## Details

The mixture residual is
`r_s = sum_{i in s} sum_k p_ik (y_i - x_i'beta_k - phi_{k,s}) / n_s`,
the responsibility-weighted residual aggregated to unit `s` (per-regime
residuals are ambiguous off-support; documented in the spec).

## References

spdep permutation Moran's I (Bivand et al.); Cliff & Ord (1981).

## Examples

``` r
# \donttest{
set.seed(1)
d <- sim_spmixqr(n = 200, G = 1, tau = 0.5, spatial_error = TRUE,
                 lattice = 6, seed = 1)
fit <- spmixqr(y ~ x, d$data, coords = d$region, G = 1, tau = 0.5,
               spatial_error = TRUE, spatial_coef = FALSE, spatial_W = d$spatial_W,
               variance = "none", control = spmixqr_control(nstart = 1L, seed = 1))
#> spatial_error = TRUE: turning the spatial gate off (the free spatial gate aliases with the CAR intercept surface). Set gating = ~z for a covariate gate.
moran_resid(fit, nsim = 199)
#> Permutation Moran's I (responsibility-weighted unit residuals)
#>   before CAR: I = 0.3320   p = 0.0300   (units = 36)
#>   after  CAR: I = -0.0908   p = 0.4450   (units = 36)
# }
```

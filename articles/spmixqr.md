# Getting started with spmixqr

`spmixqr` fits **quantile regression for spatial data while accounting
for spatial autocorrelation**. A conditional-autoregressive (CAR)
spatial random effect, built from a weights matrix, absorbs residual
spatial dependence; a Moran’s I diagnostic checks whether any spatial
structure is left over. It extends to finite mixtures of spatial
quantile regressions.

``` r
library(spmixqr)
```

## A first fit

We model the median of Columbus neighbourhood crime on income and
housing value, with a CAR spatial-error term on the queen-contiguity
graph.

``` r
data(columbus_crime); data(columbus_W)
W <- spq_weights(columbus_W, type = "supplied")

fit <- spmixqr(crime ~ income + hoval, data = columbus_crime,
               coords = factor(columbus_crime$id), G = 1, tau = 0.5,
               spatial_error = TRUE, spatial_coef = FALSE, spatial_W = W,
               control = spmixqr_control(nstart = 1L, seed = 1))
#> spatial_error = TRUE: turning the spatial gate off (the free spatial gate aliases with the CAR intercept surface). Set gating = ~z for a covariate gate.
coef(fit)
#>                regime1
#> (Intercept) 69.1571542
#> income      -2.0459672
#> hoval       -0.1137279
```

Higher neighbourhood income is associated with substantially less crime
at the median.

## Did the spatial term absorb the autocorrelation?

``` r
moran_resid(fit)
#> Permutation Moran's I (responsibility-weighted unit residuals)
#>   before CAR: I = 0.1710   p = 0.0490   (units = 49)
#>   after  CAR: I = 0.1448   p = 0.0990   (units = 49)
```

Residual spatial clustering that is present before the CAR term drops
toward non-significance after it: the spatial effect soaks up the
neighbourhood spillover.

See
[`vignette("spmixqr-primer")`](https://kvenkita.github.io/spmixqr/articles/spmixqr-primer.md)
for the full worked example (how income’s effect varies across the crime
distribution, mapping the spatial effect, inference, and diagnostics)
and `vignette("spmixqr-mixtures")` for the finite-mixture features.

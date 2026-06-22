# spmixqr ![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)

**Spatial quantile regression with autocorrelation.** `spmixqr` fits
quantile regression for spatial data and accounts for the spatial
dependence that ordinary quantile regression ignores. A
conditional-autoregressive (**CAR/ICAR**) spatial random effect — built
from a **weights matrix** you derive from queen/rook contiguity,
distance bands, k-nearest neighbours, or supply yourself — absorbs
residual spatial autocorrelation, and a permutation **Moran’s I**
diagnostic tells you whether any spatial structure is left unmodelled.
It answers questions a mean regression cannot: how the *tails* of an
outcome respond to covariates across a map, once spatial spillover is
taken into account.

With one region it is a maintained **spatial-error quantile regression**
(the lineage of Reich, Fuentes & Dunson 2011); it extends to **finite
mixtures** of spatial quantile regressions with spatially varying mixing
and covariate-effect surfaces. It is the spatial member of the `mixqr`
family and reuses that estimation core.

## Why

No maintained R package fits spatial quantile regression — the only
dedicated one, `McSpatial`, has been CRAN-archived since 2021, and
`BSquare` since 2019 — and none combines a quantile likelihood with a
CAR spatial-error term and contiguity/distance weights tooling.
`spmixqr` fills that gap.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("kvenkita/mixqr")     # dependency
remotes::install_github("kvenkita/spmixqr")
```

## Quick start

``` r
library(spmixqr)

data(columbus_crime); data(columbus_W)
W <- spq_weights(columbus_W, type = "supplied")   # or spq_weights(polys, "queen")

fit <- spmixqr(crime ~ income + hoval, data = columbus_crime,
               coords = factor(columbus_crime$id), G = 1, tau = 0.5,
               spatial_error = TRUE, spatial_coef = FALSE, spatial_W = W)
summary(fit)              # quantile coefficients + the spatial-error term
moran_resid(fit)          # residual spatial autocorrelation, before vs after the CAR term
phi_surface(fit)          # the fitted spatial-effect surface, per region
```

A fully worked walk-through (Columbus neighbourhood crime: how income’s
effect on the crime tails varies, and absorbing neighbourhood spillover)
is in
[`vignette("spmixqr-primer")`](https://kvenkita.github.io/spmixqr/articles/spmixqr-primer.md).

## Key features

- **Spatial weights, your way.**
  [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  builds `W` from **queen/rook contiguity** (`sf`/`sp` polygons),
  **distance bands**, or **k-nearest neighbours**, or accepts a
  **user-supplied** matrix.
- **CAR spatial-error term** (`spatial_error = TRUE`): a mean-zero
  conditional- autoregressive surface with a Leroux proper-CAR / ICAR
  precision that absorbs residual spatial autocorrelation, with a
  before/after **Moran’s I** diagnostic and a
  [`phi_surface()`](https://kvenkita.github.io/spmixqr/reference/phi_surface.md)
  map.
- **Quantile regression at any `tau`**, so covariate effects can differ
  across the distribution (tails vs centre) — single-region or as a
  finite **mixture** of spatial quantile regressions (spatially varying
  mixing and slope surfaces).
- Penalised EM (sparse CAR penalty; `mgcv` smooths for the gate/slope
  surfaces); exact reductions to `quantreg`, `mixqr`, and `mixqrgate`.
- Classification-aware inference: a spatial-block bootstrap
  (recommended) plus a fast sandwich;
  [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md)
  for regimes and smoothing; a full S3 surface, surface accessors, and
  maps.

## Citation

``` r
citation("spmixqr")
```

The method builds on Reich, Fuentes & Dunson (2011, *JASA*), Wu & Yao
(2016, *CSDA*), and Fernandes, Guerre & Horta (2021, *JBES*).

## License

MIT © Kailas Venkitasubramanian.

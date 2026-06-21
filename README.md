# spmixqr <img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg" align="right" />

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/kvenkita/spmixqr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kvenkita/spmixqr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Spatial finite mixtures of quantile regressions.** `spmixqr` partitions spatial
data into a few latent *regimes*, each a quantile regression, where the probability
of belonging to a regime varies smoothly over space and each regime's covariate
slopes are spatial surfaces. It answers questions a mean regression cannot: how the
*tails* of an outcome respond to covariates, and how that response, and the latent
regime structure itself, changes across a map. With one regime it is a penalised
spatially-varying-coefficient quantile regression in the lineage of Reich, Fuentes &
Dunson (2011).

It is the spatial member of the `mixqr` family ([`mixqr`](https://github.com/kvenkita/mixqr)
= mixtures of quantile regressions; [`mixqrgate`](https://github.com/kvenkita/mixqrgate)
= covariate-dependent gating), and reuses that estimation core.

## Why

No maintained R package fits spatial quantile regression (`BSquare` and `McSpatial`
are CRAN-archived), and none fits a *mixture* of quantile regressions with a spatial
gate. `spmixqr` fills both gaps.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kvenkita/mixqr")     # dependency
remotes::install_github("kvenkita/spmixqr")
```

## Quick start

```r
library(spmixqr)

d <- sim_spmixqr(n = 250, G = 2, tau = 0.5, seed = 1)
fit <- spmixqr(y ~ x, data = d$data, coords = d$coords, G = 2, tau = 0.5,
               control = spmixqr_control(nstart = 4L, seed = 1))
fit
summary(fit)

predict(fit, type = "prob")        # gate probabilities over space
coef_surface(fit)                  # spatially varying slope surface
plot(fit, which = "gate")          # map where each regime dominates
```

A real-data walk-through (soil contamination on the Meuse flood plain) is in
`vignette("spmixqr-primer")`.

## Key features

* Spatially varying mixing (softmax over a `mgcv` thin-plate / GP / MRF basis) **and**
  spatially varying component slopes.
* CAR spatial-error module (`spatial_error = TRUE`): a per-regime mean-zero
  conditional-autoregressive surface `phi_k` with a Leroux proper-CAR/ICAR precision,
  built from a `spq_weights()` matrix (queen/rook contiguity, distance, k-NN, or
  user-supplied). Includes a before/after residual Moran's I diagnostic.
* Penalised EM reusing `mixqr`; convolution-smoothed component step so a roughness
  penalty applies; exact reduction to weighted quantile regression for flat slopes.
* Classification-aware inference: spatial-block / xy bootstrap (recommended) plus a
  fast classification-conditional sandwich.
* `spmixqr_select()` for the number of regimes and smoothing; full S3 surface;
  surface accessors and maps; calibrated diagnostics including a label-stability
  warning.
* Reductions to `quantreg`, `mixqr`, and `mixqrgate` are tested exactly.

## Citation

```r
citation("spmixqr")
```

The method builds on Reich, Fuentes & Dunson (2011, *JASA*), Wu & Yao (2016, *CSDA*),
and Fernandes, Guerre & Horta (2021, *JBES*).

## License

MIT © Kailas Venkitasubramanian.

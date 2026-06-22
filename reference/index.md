# Package index

## Model fitting

Fit a spatial quantile regression (single-region or mixture) and
configure the EM.

- [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  : Fit a spatial finite mixture of quantile regressions

- [`spmixqr_control()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_control.md)
  :

  Control parameters for
  [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)

- [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md)
  : Select regimes and smoothing parameters

## Spatial weights and autocorrelation (CAR / NNGP error)

Build the weights matrix or an NNGP Matern precision (point data),
diagnose residual spatial autocorrelation, and read the spatial effect.

- [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  : Construct a spatial weights object for the CAR-error term
- [`moran_resid()`](https://kvenkita.github.io/spmixqr/reference/moran_resid.md)
  : Permutation Moran's I on the spatial-unit residuals of an spmixqr
  fit
- [`phi_surface()`](https://kvenkita.github.io/spmixqr/reference/phi_surface.md)
  : CAR spatial-error surfaces (with uncertainty)

## Mixture surfaces and prediction

- [`spmixqr_basis()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_basis.md)
  : Build a spatial basis and its roughness penalty
- [`coef_surface()`](https://kvenkita.github.io/spmixqr/reference/coef_surface.md)
  : Spatial coefficient surfaces
- [`gate_surface()`](https://kvenkita.github.io/spmixqr/reference/gate_surface.md)
  : Spatial gate surfaces
- [`predict(`*`<spmixqr>`*`)`](https://kvenkita.github.io/spmixqr/reference/predict.spmixqr.md)
  : Predictions from a spatial quantile mixture model
- [`plot(`*`<spmixqr>`*`)`](https://kvenkita.github.io/spmixqr/reference/plot.spmixqr.md)
  : Plot a spatial quantile mixture fit

## Simulation and data

- [`sim_spmixqr()`](https://kvenkita.github.io/spmixqr/reference/sim_spmixqr.md)
  : Simulate from a spatial quantile mixture model
- [`sim_spatial_confound()`](https://kvenkita.github.io/spmixqr/reference/sim_spatial_confound.md)
  : Simulate point data with spatial confounding (for Spatial+ / NNGP
  demos)
- [`columbus_crime`](https://kvenkita.github.io/spmixqr/reference/columbus_crime.md)
  : Columbus, Ohio neighbourhood crime
- [`columbus_W`](https://kvenkita.github.io/spmixqr/reference/columbus_W.md)
  : Queen-contiguity weights for the Columbus neighbourhoods
- [`lucas_house`](https://kvenkita.github.io/spmixqr/reference/lucas_house.md)
  : Lucas County (Ohio) house sales (subsample)
- [`nc_sids`](https://kvenkita.github.io/spmixqr/reference/nc_sids.md) :
  North Carolina SIDS counts (areal census example)
- [`nc_sids_W`](https://kvenkita.github.io/spmixqr/reference/nc_sids_W.md)
  : Queen-contiguity weights for the North Carolina SIDS counties
- [`meuse_zinc`](https://kvenkita.github.io/spmixqr/reference/meuse_zinc.md)
  : Meuse river zinc concentrations

## Package

- [`spmixqr-package`](https://kvenkita.github.io/spmixqr/reference/spmixqr-package.md)
  : spmixqr: Spatial Finite Mixtures of Quantile Regressions

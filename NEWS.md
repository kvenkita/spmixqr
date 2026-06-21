# spmixqr 0.2.4

More interpretable spatial-error (CAR) coefficients.

* `phi_surface()` gains `ci = TRUE` (adds bootstrap `se`, `lower`, `upper`, and a
  `credible` flag for units whose interval excludes zero) and `scale = "exp"` (adds
  `mult = exp(phi)`, the multiplicative deviation for log-outcome models, e.g. 1.18 =
  "about 18% above expected"). Standard errors propagate the fit's coefficient
  covariance through the sum-to-zero constraint, `Var(phi) = T V T'`; the bootstrap is
  recommended over the sandwich, which ignores the penalty.
* `plot(fit, which = "phi", credible = TRUE)` greys out units that are not reliably
  different from zero, mapping only the credible hot/cold spots.
* The primer now reads the spatial effect with uncertainty (a credible-hotspot map and
  an outcome-scale interpretation); the spatial-error article shows the `exp()`
  multiplicative reading on log zinc.

# spmixqr 0.2.3

Mixture example moved to a larger, appropriate dataset. Documentation/data only.

* New built-in dataset `lucas_house`: a 2000-row subsample of Lucas County, Ohio house
  sales. The primer's mixture example (Example B) now fits a two-regime spatially gated
  mixture of quantile regressions there (housing submarkets that discount age at
  different rates), with point maps of the gate-probability surface, the assigned
  submarket, and the residuals. The 49-neighbourhood Columbus data is too small to
  resolve a mixture and is no longer used for it (it remains Example A, the single-region
  CAR fit).

# spmixqr 0.2.2

Primer expanded with maps and a mixture example. Documentation only.

* The Columbus shapefile now ships in `inst/shapes/columbus.gpkg`. The primer reads it
  with `sf::st_read()` and draws **choropleth maps** of the spatial-error surface and the
  residuals (replacing the earlier centroid scatter), demonstrating the
  shapefile-to-`sf`-to-map workflow.
* A second worked example fits a **two-regime mixture** (`G = 2`) of quantile
  regressions on Columbus crime, maps the assigned regimes and residuals, and reports
  the per-regime coefficient table.

# spmixqr 0.2.1

Documentation refocus: the spatial-dependence features (CAR spatial-error term,
weights-matrix construction, Moran's I autocorrelation diagnostics) are now the
package's primary framing, with the finite-mixture features presented as an
extension. No change to the estimation code or the API.

* New built-in dataset `columbus_crime` (49 Columbus, Ohio neighbourhoods) and its
  queen-contiguity weights `columbus_W`, for the rewritten primer.
* The primer (`vignette("spmixqr-primer")`) is rewritten as a fully worked
  spatial-error quantile regression on Columbus crime: building the weights, fitting,
  reading the income effect across the crime distribution, absorbing residual spatial
  autocorrelation (Moran's I before/after), mapping the spatial effect, and reporting.
* The earlier mixture/gate/slope walk-through moves to `vignette("spmixqr-mixtures")`.
* README, DESCRIPTION, and the pkgdown site lead with the spatial-autocorrelation
  capabilities.

# spmixqr 0.2.0

Adds a conditional-autoregressive (CAR) spatial-**error** module. Fully
back-compatible: `spatial_error = FALSE` (the default) reproduces v0.1.0 exactly.

* New `spatial_error = TRUE` argument to `spmixqr()`: each regime carries a mean-zero
  per-unit CAR/ICAR spatial random effect `phi_k` with a Leroux proper precision
  `Q(alpha) = alpha (D - W) + (1 - alpha) I` (or intrinsic `D - W`). A per-connected-
  component sum-to-zero constraint is absorbed into the incidence/precision (à la
  `mgcv`'s `absorb.cons`), so the intercept carries the level and `phi_k` is identified.
  At `G = 1` this is a maintained single-population spatial-error quantile regression;
  at `G > 1` (covariate/constant gate) it is a mixture of spatial quantile regressions.
* New weights layer `spq_weights()`: builds a symmetric, sparse `W` from `sf`/`sp`
  polygons (queen/rook contiguity), point coordinates (distance band, k-NN, symmetrised),
  an `spdep` `nb`/`listw`, or a user matrix (`supplied`). `make_car_precision()` forms
  the CAR/ICAR precision; `print.spq_weights()`.
* The component M-step is rebranched: when `spatial_error = TRUE` it routes through a
  new **sparse** penalised-Newton solver (each Newton step solved with `Matrix::solve()`
  on the sparse `(p + L')` system) on the augmented `[X, R]` design, regardless of
  `spatial_coef`, so the CAR penalty is never dropped. The dense solver remains the
  non-CAR path.
* Identification guardrail: `spatial_error = TRUE` turns the spatial gate off by default
  (with a message); an explicit `spatial_gate = TRUE` alongside `spatial_error = TRUE`
  raises an explanatory error (a free spatial gate aliases with the CAR level surface).
* `moran_resid()`: permutation Moran's I on the responsibility-weighted residual
  aggregated to the spatial unit, reported before/after the CAR term. `phi_surface()`
  tidy accessor; `plot(which = "phi")`; a spatial-error block in `summary()`.
* `alpha` is fixed (default proper-CAR `0.95`, disclosed as weakly identified, not
  selected by the check loss); `lambda_error` (= `lambda_phi`) is selected by BIC with
  a CAR effective-df term. `spmixqr_select()` gains `lambda_error`/`car`/`car_alpha`
  pass-through.
* Inference: `bootstrap_vcov()` gains an areal connected-component block mode
  (`spdep::n.comp.nb`) for CAR fits (the point path silently fell back to i.i.d.).
* `sim_spmixqr(spatial_error = TRUE)` generates lattice data with a known CAR effect.
* New shipped areal example `nc_sids` (North Carolina SIDS, 100 counties) with queen
  contiguity weights `nc_sids_W`, built offline (CSV/`.rda`-decoupled) from `sf`.
* Imports `spdep`, `Matrix`, `methods` (the sparse CAR precision/solver); `sf` in
  Suggests (polygon input).

# spmixqr 0.1.0

First release.

* `spmixqr()` fits spatial finite mixtures of quantile regressions: latent regimes
  with spatially varying mixing probabilities (a softmax over a spatial basis) and
  quantile-regression components whose covariate slopes are spatial surfaces (scalar
  per-regime intercepts). With `G = 1` it is a penalised spatially-varying-coefficient
  quantile regression in the lineage of Reich, Fuentes & Dunson (2011).
* Penalised EM reusing the `mixqr` component machinery (`weighted_rq`,
  `constrained_kde`) and a vendored, spatially penalised multinomial-logit gate. The
  spatially-varying-slope component step uses a convolution-smoothed check loss
  (Fernandes, Guerre & Horta 2021) with bandwidth annealing, so a roughness penalty
  applies; it reduces to exact weighted quantile regression when slopes are flat.
* Spatial bases (thin-plate, low-rank Gaussian process, Markov random field) and
  roughness penalties built with `mgcv`; coordinates standardised internally and the
  transform stored for prediction at new locations.
* Classification-aware inference: a fast classification-conditional sandwich, and a
  recommended spatial-block / xy-pairs bootstrap with cross-replicate label alignment.
* Identifiability and labels: a free spatial intercept surface is not separately
  identified from the gate, so spatial level is carried by the gate; components are
  ordered by a sign-stable separating slope, with a `label_stability` diagnostic that
  warns when slope surfaces cross.
* `spmixqr_select()` (BIC or cross-validated check loss), full S3 method surface,
  `coef_surface()` / `gate_surface()` accessors, `sim_spmixqr()` generator, and the
  built-in `meuse_zinc` dataset.
* Validation: recovery at held-out locations, Monte-Carlo coverage, two negative
  controls (all-flat and aliasing-targeted), and reductions to `quantreg`, `mixqr`,
  and `mixqrgate` (see `inst/benchmarks` and `inst/replication`).

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

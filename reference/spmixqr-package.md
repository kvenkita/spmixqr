# spmixqr: Spatial Finite Mixtures of Quantile Regressions

Fits spatial finite mixtures of quantile regressions. The mixing
(gating) probabilities vary over space through a low-rank spatial basis,
and each regime's covariate-effect slopes are spatially varying surfaces
(with a scalar per-regime intercept). Estimation is a penalised
expectation-maximisation algorithm that reuses the mixqr component
machinery and a vendored, spatially penalised multinomial-logit gate.
The single-regime case is a penalised spatially-varying-coefficient
quantile regression in the lineage of Reich, Fuentes and Dunson (2011).

## Entry points

- [`spmixqr()`](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  — fit a spatial mixture of quantile regressions.

- [`spmixqr_select()`](https://kvenkita.github.io/spmixqr/reference/spmixqr_select.md)
  — choose the number of regimes and/or the smoothing parameters.

- [`sim_spmixqr()`](https://kvenkita.github.io/spmixqr/reference/sim_spmixqr.md)
  — simulate from the model (for validation).

- [`coef_surface()`](https://kvenkita.github.io/spmixqr/reference/coef_surface.md),
  [`gate_surface()`](https://kvenkita.github.io/spmixqr/reference/gate_surface.md),
  [`phi_surface()`](https://kvenkita.github.io/spmixqr/reference/phi_surface.md)
  — spatial-surface accessors (slope, gate, and CAR spatial-error
  surfaces).

- [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  — build the CAR spatial weights matrix.

- [`moran_resid()`](https://kvenkita.github.io/spmixqr/reference/moran_resid.md)
  — residual Moran's I before/after the CAR term.

## References

Reich, B. J., Fuentes, M. and Dunson, D. B. (2011). Bayesian spatial
quantile regression. *Journal of the American Statistical Association*
106, 6–20.

Wu, C. and Yao, W. (2016). Mixtures of quantile regressions.
*Computational Statistics & Data Analysis* 93, 162–176.

Fernandes, M., Guerre, E. and Horta, E. (2021). Smoothing quantile
regressions. *Journal of Business & Economic Statistics* 39, 338–357.

## See also

Useful links:

- <https://github.com/kvenkita/spmixqr>

- <https://kvenkita.github.io/spmixqr/>

- Report bugs at <https://github.com/kvenkita/spmixqr/issues>

## Author

**Maintainer**: Kailas Venkitasubramanian <kailasv@gmail.com>

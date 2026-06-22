# Spatial label-coherence diagnostic.

Given each regime's ordering-covariate slope evaluated at every location
(`slope_surface`, an `n x G` matrix), returns the fraction of locations
where the pointwise ranking of the regimes disagrees with the global
(mean-slope) ranking. A high value means slope surfaces cross in space
and a single global label is not pointwise coherent (v1 limitation;
per-region alignment is v2).

## Usage

``` r
label_stability(slope_surface)
```

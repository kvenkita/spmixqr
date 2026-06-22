# Bootstrap covariance (xy-pairs or spatial-block) of the EM pipeline.

Reuses the fixed training basis rows (resampling observations, not
rebuilding the basis), refits, aligns each replicate's regimes to the
point estimate by matching constant coefficients, and returns the
empirical covariance of the stacked gate and coefficient parameters.
Under-covers with few blocks (Lahiri 2003); a minimum-sites guard
applies.

## Usage

``` r
bootstrap_vcov(obj, data = NULL, B = NULL, block = NULL)
```

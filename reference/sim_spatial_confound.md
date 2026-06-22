# Simulate point data with spatial confounding (for Spatial+ / NNGP demos)

Generates point-referenced data in which a covariate is correlated with
a smooth spatial field that also drives the outcome – the
spatial-confounding setting of Dupont, Wood & Augustin (2022). A
high-frequency field `f(s)` is used so the covariate's spatial signal
out-resolves a penalised spatial-error term (the resolution gap that
makes Spatial+ effective; Frisch–Waugh–Lovell). With `confound = 0` it
is a negative control (no confounding). The error is symmetric, so the
true `tau`-quantile slope equals `beta` at every `tau`.

## Usage

``` r
sim_spatial_confound(
  n = 400L,
  beta = 1,
  confound = 0.9,
  amp = 1.2,
  freq = 6,
  sd = 0.5,
  seed = NULL
)
```

## Arguments

- n:

  Number of points.

- beta:

  True covariate slope (recovered by a deconfounded estimator).

- confound:

  Strength of covariate–space correlation (`0` = negative control).

- amp:

  Amplitude of the spatial field in the outcome.

- freq:

  Spatial frequency of the field (higher = finer; default 6).

- sd:

  Error scale.

- seed:

  Optional seed.

## Value

A list with `data` (`y`, `x`), `coords` (`n x 2`), and `truth` (`beta`,
the field `f`, `confound`).

## References

Dupont, Wood & Augustin (2022, Biometrics).

## Examples

``` r
d <- sim_spatial_confound(n = 200, seed = 1)
str(d$truth)
#> List of 3
#>  $ beta    : num 1
#>  $ f       : num [1:200] 0.966 1.045 -1.29 -0.783 1.401 ...
#>  $ confound: num 0.9
```

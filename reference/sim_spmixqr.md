# Simulate from a spatial quantile mixture model

Generates point (or areal) data from a two-regime spatial mixture of
quantile regressions: a smooth spatial gate and spatially varying
component slopes (à la Reich, Fuentes & Dunson 2011, eq. 17, inside a
mixture). Used for the validation study. The conditional `tau`-quantile
is exact by construction (the error's `tau`-quantile is subtracted).

## Usage

``` r
sim_spmixqr(
  n = 200L,
  G = 2L,
  tau = 0.5,
  gate_slope = 2.5,
  coef_slope = 1.5,
  sep = 3,
  sd = 0.7,
  crossing = FALSE,
  error = c("normal", "ald"),
  spatial_error = FALSE,
  lattice = 8L,
  car_alpha = 0.95,
  car_rho = 1.5,
  seed = NULL
)
```

## Arguments

- n:

  Number of observations / distinct locations.

- G:

  Number of regimes (2 or 3 supported).

- tau:

  Quantile level.

- gate_slope:

  Strength of the spatial gate (0 = constant gate; negative control uses
  0).

- coef_slope:

  Strength of the spatial slope variation (0 = flat slopes).

- sep:

  Regime separation (intercept gap).

- sd:

  Error scale.

- crossing:

  If `TRUE`, make slope surfaces cross in space (stresses labels).

- error:

  `"normal"` or `"ald"` component error.

- spatial_error:

  If `TRUE`, generate data on a regular `lattice x lattice` grid and add
  a per-regime CAR spatial random effect `phi_k` drawn from the
  proper-CAR precision built from the lattice's rook contiguity. The
  returned list then gains `spatial_W` (the
  [`spq_weights()`](https://kvenkita.github.io/spmixqr/reference/spq_weights.md)
  object), and `truth$phi` (the L x G true CAR surface, each column
  mean-zero) and `truth$unit` (observation -\> unit).

- lattice:

  Side length of the square lattice when `spatial_error = TRUE`
  (`L = lattice^2` units); `n` observations are assigned to units (with
  repeats).

- car_alpha, car_rho:

  Proper-CAR alpha used to build the precision and the marginal SD
  scaling of the simulated CAR effect (`car_rho` = phi SD).

- seed:

  Optional seed.

## Value

A list with `data` (y, x, x2), `coords` (n x 2), and `truth`. With
`spatial_error = TRUE` also `spatial_W` and a `region` column /
`truth$phi`.

## Examples

``` r
d <- sim_spmixqr(n = 150, seed = 1)
str(d$truth)
#> List of 9
#>  $ z         : int [1:150] 2 2 2 2 1 2 2 2 2 1 ...
#>  $ pi        : num [1:150, 1:2] 0.813 0.719 0.484 0.149 0.857 ...
#>  $ q         : num [1:150] 2.367 2.524 4.895 5.933 -0.367 ...
#>  $ slope_x   : num [1:150] -2.398 -2.558 -2.859 -3.362 0.803 ...
#>  $ intercepts: num [1:2] 0 3
#>  $ tau       : num 0.5
#>  $ G         : int 2
#>  $ gate_slope: num 2.5
#>  $ coef_slope: num 1.5
```

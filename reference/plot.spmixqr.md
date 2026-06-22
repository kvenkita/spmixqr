# Plot a spatial quantile mixture fit

Base-graphics maps and diagnostics. Richer ggplot2 maps are shown in the
package primer.

## Usage

``` r
# S3 method for class 'spmixqr'
plot(
  x,
  which = c("gate", "coef", "phi", "class", "density"),
  regime = NULL,
  credible = FALSE,
  level = 0.95,
  ...
)
```

## Arguments

- x:

  A fitted
  [spmixqr](https://kvenkita.github.io/spmixqr/reference/spmixqr.md)
  object.

- which:

  One of `"gate"` (mixing-probability map), `"coef"` (slope-surface
  map), `"phi"` (CAR spatial-error surface), `"class"` (classified map),
  or `"density"` (component error densities).

- regime:

  Which regime to map for `"gate"`/`"coef"`/`"phi"` (default 2 for gate,
  1 for coef/phi).

- credible:

  For `which = "phi"`, show only the **credible** hot/cold spots (units
  whose `level` confidence interval for the CAR effect excludes zero);
  other units are greyed. Requires a fit carrying a covariance
  (`variance = "boot"` or `"sandwich"`).

- level:

  Confidence level used when `credible = TRUE`.

- ...:

  Passed to the underlying plotting call.

## Value

Invisibly `NULL`.

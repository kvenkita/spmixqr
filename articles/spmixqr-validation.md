# Validation: recovery, coverage, and negative controls

``` r
library(spmixqr)
set.seed(1)
```

This article summarises the validation evidence. A fuller grid is
shipped in `inst/benchmarks/coverage.R` and
`inst/replication/sim_grid.R`.

## Recovery at held-out locations

We simulate a two-regime spatial mixture, fit on a training set, and
measure the recovered gate surface against the truth at held-out
locations.

``` r
d <- sim_spmixqr(n = 400, G = 2, tau = 0.5, seed = 2)
tr <- seq_len(300); te <- 301:400
fit <- spmixqr(y ~ x, d$data[tr, ], coords = d$coords[tr, ], G = 2, tau = 0.5,
               variance = "none", control = spmixqr_control(nstart = 4L, seed = 1))
pr_te <- predict(fit, newdata = d$data[te, ], newcoords = d$coords[te, ], type = "prob")
gate_cor <- cor(pr_te[, 2], d$truth$pi[te, 2])
cl <- apply(predict(fit, newdata = d$data[te, ], newcoords = d$coords[te, ],
                    type = "prob"), 1, which.max)
acc <- max(mean(cl == d$truth$z[te]), mean((3 - cl) == d$truth$z[te]))
round(c(held_out_gate_cor = abs(gate_cor), held_out_accuracy = acc), 3)
#> held_out_gate_cor held_out_accuracy 
#>             0.994             0.750
```

The gate surface and the regime membership are recovered at locations
not used in fitting.

## Negative control

With no spatial structure, the method must not invent any. We simulate a
constant gate and flat slopes, fit with a strong penalty, and confirm
the gate surface is essentially flat.

``` r
nc <- sim_spmixqr(n = 250, gate_slope = 0, coef_slope = 0, seed = 3)
fit_nc <- spmixqr(y ~ x, nc$data, coords = nc$coords, G = 2, tau = 0.5,
                  lambda_gate = 10, lambda_coef = 10, variance = "none",
                  control = spmixqr_control(nstart = 4L, seed = 1))
gs <- gate_surface(fit_nc)
round(sd(gs$prob[gs$regime == "2"]), 3)   # near 0 => no manufactured structure
#> [1] 0.106
```

## Inference

Coverage is studied in `inst/benchmarks/coverage.R`. The recommended
bootstrap attains near-nominal coverage in well-spanned designs; the
classification-conditional sandwich is faster but optimistic, and is
disclosed as such (Reich, Fuentes, and Dunson 2011). The
asymmetric-Laplace component likelihood is a working likelihood, so we
default to the bootstrap for reporting.

## References

Reich, Brian J., Montserrat Fuentes, and David B. Dunson. 2011.
“Bayesian Spatial Quantile Regression.” *Journal of the American
Statistical Association* 106 (493): 6–20.
<https://doi.org/10.1198/jasa.2010.ap09237>.

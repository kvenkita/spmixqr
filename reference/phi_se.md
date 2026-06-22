# Standard errors of the CAR spatial effect phi (L x G).

Propagates the fit's stored coefficient covariance through the
sum-to-zero constraint transform: `Var(phi_k) = T V_red,k T'`, with
`V_red,k` the CAR sub-block of the regime-k coefficient covariance and
`T` the constraint-absorption basis. The **bootstrap**
(`variance = "boot"`) is recommended: it refits the penalised pipeline
and so reflects the shrinkage of the random effect, whereas the
classification- conditional sandwich is a fast alternative that ignores
the penalty and can disagree. With few spatial blocks the bootstrap
intervals can be optimistic (small-sample caveat). Returns `NULL` if no
covariance is stored.

## Usage

``` r
phi_se(object)
```

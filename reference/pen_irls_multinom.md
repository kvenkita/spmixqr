# Fit a penalised weighted multinomial logit by Newton/IRLS.

Fit a penalised weighted multinomial logit by Newton/IRLS.

## Usage

``` r
pen_irls_multinom(Z, P, Pen, maxit = 50L, tol = 1e-08)
```

## Arguments

- Z:

  Gate design `n x q1`.

- P:

  Fractional responses `n x G` (responsibilities; rows sum to 1).

- Pen:

  Per-class penalty matrix `q1 x q1` (PSD).

- maxit, tol:

  Newton iterations / tolerance.

## Value

list(gamma `q1 x (G-1)`, hessian, pi, converged).

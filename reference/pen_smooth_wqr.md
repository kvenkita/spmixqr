# Penalised smoothed weighted quantile regression (one regime's component step).

Solves the penalised smoothed objective by penalised Newton with
**bandwidth annealing** (a wide bandwidth is reduced geometrically to
the target, warm- starting each stage) and step-halving. Annealing is
essential: at a tiny target bandwidth the kernel Hessian is
near-singular and undamped Newton diverges (verified). At convergence
the solution approximates the target-bandwidth minimiser, and as the
target bandwidth shrinks it approaches exact weighted quantile
regression.

## Usage

``` r
pen_smooth_wqr(
  Xt,
  y,
  tau,
  w,
  Pen,
  h,
  floor_abs = 0.001,
  beta_init = NULL,
  maxit = 50L,
  tol = 1e-07
)
```

## Arguments

- Xt:

  Augmented component design `n x P` (`[1, x_j, x_j*B(s)]`, j\>=2).

- y:

  Response.

- tau:

  Quantile level.

- w:

  Observation weights (responsibilities) length n.

- Pen:

  `P x P` penalty matrix (PSD).

- h:

  Bandwidth *rate* (dimensionless); the absolute bandwidth is this rate
  times the residual scale, floored at `floor_abs`.

- floor_abs:

  Absolute minimum bandwidth.

- beta_init:

  Optional warm start (ignored if `NULL` or all-zero).

- maxit, tol:

  Newton iterations / tolerance per annealing stage.

## Value

list(beta, hessian, fitted, h (absolute bandwidth used), converged).

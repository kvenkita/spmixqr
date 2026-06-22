# Sparse penalised smoothed weighted quantile regression for the CAR design.

Identical annealed penalised-Newton logic to
[`pen_smooth_wqr()`](https://kvenkita.github.io/spmixqr/reference/pen_smooth_wqr.md),
but the augmented design `[X, Rt]` (with `Rt` the constraint-absorbed
CAR incidence) and the penalty block are stored as Matrix sparse
objects, and the `(p + L')^2` Newton system is assembled and solved
sparsely via
[`Matrix::solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html) /
a sparse Cholesky. The dense
[`pen_smooth_wqr()`](https://kvenkita.github.io/spmixqr/reference/pen_smooth_wqr.md)
stays the non-CAR path; this solver is the CAR path
(`spatial_error = TRUE`). A dense fallback is used for tiny `L'`.

## Usage

``` r
pen_smooth_wqr_car(
  Xt,
  y,
  tau,
  w,
  Pen,
  h,
  floor_abs = 0.001,
  beta_init = NULL,
  maxit = 50L,
  tol = 1e-07,
  dense_max = 40L
)
```

## Arguments

- Xt:

  Augmented component design, a sparse `Matrix` `n x P` (`P = p + L'`).

- y:

  Response.

- tau:

  Quantile level.

- w:

  Observation weights (responsibilities) length n.

- Pen:

  `P x P` sparse penalty matrix (PSD; zero on the beta block).

- h:

  Bandwidth rate.

- floor_abs:

  Absolute minimum bandwidth.

- beta_init:

  Optional warm start.

- maxit, tol:

  Newton iterations / tolerance per annealing stage.

- dense_max:

  If `P <= dense_max`, fall back to the dense solver path.

## Value

list(beta, hessian, fitted, h, converged) (beta a plain numeric vector).

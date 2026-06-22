# Constraint-absorbed incidence and precision for the CAR block.

Imposes a per-connected-component sum-to-zero constraint on `phi` (so
`beta`'s intercept carries the level and `phi` is the mean-zero spatial
deviation), absorbing the constraint into a reduced incidence basis `Rt`
(`n x L'`, `L' = L - n_comp`) and reduced precision `Qt` (`L' x L'`),
analogous to mgcv's `absorb.cons`. A tiny epsilon-ridge is added to `Qt`
for conditioning. The CAR effect on the original `L` units is recovered
as `phi = Tmat %*% phi_reduced`, where `Tmat` (`L x L'`) is the
constraint null-space basis stored in the return value.

## Usage

``` r
absorb_car_constraint(R, Q, membership, eps = 1e-06, constrain = TRUE)
```

## Arguments

- R:

  Full unit-incidence matrix `n x L` (`R[i, s(i)] = 1`), sparse.

- Q:

  Full `L x L` CAR precision from
  [`make_car_precision()`](https://kvenkita.github.io/spmixqr/reference/make_car_precision.md).

- membership:

  Length-L connected-component labels (from
  [`components_from_W()`](https://kvenkita.github.io/spmixqr/reference/components_from_W.md)).

- eps:

  Ridge added to the reduced precision.

- constrain:

  Apply the per-component sum-to-zero (default `TRUE`, for the
  rank-deficient ICAR/CAR precision). `FALSE` for a proper full-rank
  precision (NNGP): `phi` is fit unconstrained and `Tmat` is the
  identity.

## Value

A list with `Rt` (constraint-absorbed incidence, sparse `n x L'`), `Qt`
(reduced precision, sparse `L' x L'`), `Tmat` (`L x L'` recovery basis),
and `Lp` (`L'`).

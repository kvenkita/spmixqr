# Build a CAR / ICAR precision matrix from a weights object

Forms the proper Leroux precision
`Q(alpha) = alpha (D - W) + (1 - alpha) I` (PSD; PD for `alpha < 1`), or
the intrinsic CAR precision `D - W` plus a small epsilon-ridge guard
(`car = "icar"`). The proper default needs no ridge; the ICAR precision
is rank-deficient (null dimension = number of connected components) and
the ridge only makes it numerically invertible - identification of the
level comes from the per-component sum-to-zero constraint absorbed in
[`build_designs()`](https://kvenkita.github.io/spmixqr/reference/build_designs.md).

## Usage

``` r
make_car_precision(spqw, alpha = 0.95, car = c("proper", "icar"), eps = 1e-06)
```

## Arguments

- spqw:

  An `spq_weights` object (or a list with `W`, `D`).

- alpha:

  Proper-CAR spatial-dependence strength in `[0, 1]` (default `0.95`).

- car:

  `"proper"` (Leroux) or `"icar"` (intrinsic).

- eps:

  Epsilon-ridge added to the ICAR precision for conditioning.

## Value

A symmetric sparse precision matrix (`dsCMatrix`).

## References

Leroux et al. (2000).

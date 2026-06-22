# Connected-component count and membership from a symmetric weights matrix.

Pure-`Matrix`/base graph traversal (no spdep dependency) so the
`supplied` path and the constraint-absorption code can find connected
components of any `W`.

## Usage

``` r
components_from_W(W)
```

## Arguments

- W:

  A symmetric sparse/dense weights matrix.

## Value

A list with `nc` (number of components) and `membership` (length-L
integer component label per unit).

# Resolve coordinates / areal specification.

For the CAR spatial-error path with a supplied `spatial_W` (areal data
without an spdep `nb`), `coords` is treated as a length-n vector of
region labels even when `areal` is `NULL`.

## Usage

``` r
resolve_coords(coords, areal, data, n, spatial_error = FALSE, spatial_W = NULL)
```

# Dimensionless conquer-style bandwidth rate ~ scale \* (log n / n)^(1/5). The absolute bandwidth used in the M-step is this rate times the residual scale (so the smoothing is scale-equivariant), with an absolute floor applied there.

Dimensionless conquer-style bandwidth rate ~ scale \* (log n / n)^(1/5).
The absolute bandwidth used in the M-step is this rate times the
residual scale (so the smoothing is scale-equivariant), with an absolute
floor applied there.

## Usage

``` r
smooth_bandwidth(n, scale = 1)
```

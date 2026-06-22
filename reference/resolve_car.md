# Resolve the CAR weights / precision / constraint-absorbed incidence.

Builds the `spq_weights`, the unit map (observation -\> unit), the full
incidence, the precision, and the per-component sum-to-zero constraint
absorption.

## Usage

``` r
resolve_car(spatial_W, geo, car, car_alpha, n)
```

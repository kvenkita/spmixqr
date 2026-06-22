# Lucas County (Ohio) house sales (subsample)

A seeded random subsample of 2000 single-family house sales in Lucas
County, Ohio (1993–1998), used to illustrate a finite **mixture** of
spatial quantile regressions (housing submarkets) on point data. Price
is strongly right-skewed, so the tails and latent submarkets are worth
modelling.

## Usage

``` r
lucas_house
```

## Format

A data frame with 2000 rows (sales) and 6 variables:

- price:

  Sale price (USD).

- tla:

  Total living area (square feet).

- age:

  Age of the dwelling at sale (in centuries; range about 0–1.4).

- rooms:

  Number of rooms.

- x, y:

  Projected coordinates (feet).

## Source

A subsample (seed 7) of the `house` data in the spData package (Lucas
County Auditor; see also Pace and LeSage). Built offline; see
`data-raw/make_lucas.R`.

## Examples

``` r
data(lucas_house)
summary(lucas_house$price)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    3400   41000   66000   79571   99125  650000 
```

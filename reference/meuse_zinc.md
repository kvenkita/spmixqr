# Meuse river zinc concentrations

Topsoil heavy-metal concentrations on a flood plain of the river Meuse,
near Stein (NL). A classic spatial dataset: zinc is strongly
right-skewed and its distribution shifts with distance to the river and
flood frequency, making it a natural illustration of spatially-varying
conditional quantiles and latent flood-frequency regimes.

## Usage

``` r
meuse_zinc
```

## Format

A data frame with 155 rows and 6 variables:

- zinc:

  Topsoil zinc concentration (ppm).

- dist:

  Normalised distance to the river Meuse (in `[0, 1]`).

- elev:

  Relative elevation above the channel bed (m).

- ffreq:

  Flood-frequency class (factor: 1 = once in two years, 2 = once in ten
  years, 3 = once in fifty years).

- x, y:

  Easting and northing (Rijksdriehoek / RD New, metres).

## Source

The `meuse` data from the sp package (Pebesma & Bivand), originally from
Burrough & McDonnell (1998), *Principles of Geographical Information
Systems*. Redistributed here as a built-in for reproducible examples.

## Examples

``` r
data(meuse_zinc)
summary(meuse_zinc$zinc)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   113.0   198.0   326.0   469.7   674.5  1839.0 
```

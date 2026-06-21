#' Meuse river zinc concentrations
#'
#' Topsoil heavy-metal concentrations on a flood plain of the river Meuse, near
#' Stein (NL). A classic spatial dataset: zinc is strongly right-skewed and its
#' distribution shifts with distance to the river and flood frequency, making it a
#' natural illustration of spatially-varying conditional quantiles and latent
#' flood-frequency regimes.
#'
#' @format A data frame with 155 rows and 6 variables:
#' \describe{
#'   \item{zinc}{Topsoil zinc concentration (ppm).}
#'   \item{dist}{Normalised distance to the river Meuse (in `[0, 1]`).}
#'   \item{elev}{Relative elevation above the channel bed (m).}
#'   \item{ffreq}{Flood-frequency class (factor: 1 = once in two years, 2 = once in
#'     ten years, 3 = once in fifty years).}
#'   \item{x, y}{Easting and northing (Rijksdriehoek / RD New, metres).}
#' }
#' @source The `meuse` data from the \pkg{sp} package (Pebesma & Bivand), originally
#'   from Burrough & McDonnell (1998), \emph{Principles of Geographical Information
#'   Systems}. Redistributed here as a built-in for reproducible examples.
#' @examples
#' data(meuse_zinc)
#' summary(meuse_zinc$zinc)
"meuse_zinc"

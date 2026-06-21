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

#' North Carolina SIDS counts (areal census example)
#'
#' Sudden-infant-death (SIDS) counts and live births for the 100 counties of North
#' Carolina, 1974--78, the canonical areal/lattice spatial dataset. Shipped as a tidy
#' data frame with county centroids and a derived log SIDS rate, to illustrate the CAR
#' spatial-error quantile regression on contiguous (polygon) areal units. The matching
#' queen-contiguity weights are in [nc_sids_W].
#'
#' @format A data frame with 100 rows (counties) and 9 variables:
#' \describe{
#'   \item{county}{County name (the spatial unit id; matches `nc_sids_W` row/col names).}
#'   \item{births}{Live births, 1974--78 (`BIR74`).}
#'   \item{sid}{SIDS deaths, 1974--78 (`SID74`).}
#'   \item{nwbirths}{Non-white live births, 1974--78 (`NWBIR74`).}
#'   \item{east, north}{County-centroid longitude / latitude.}
#'   \item{sids_rate}{SIDS per 1000 births (`1000 * (sid + 0.5) / (births + 1)`).}
#'   \item{log_sids}{`log(sids_rate)` -- the modelled response.}
#'   \item{pnw}{Proportion of non-white births (`nwbirths / (births + 1)`).}
#' }
#' @source The `nc.shp` North Carolina SIDS shapefile distributed with the \pkg{sf}
#'   package (originally Cressie & Read 1985; Cressie 1993, \emph{Statistics for
#'   Spatial Data}). Built offline from `sf::st_read(system.file("shape/nc.shp",
#'   package = "sf"))` and redistributed as a built-in for reproducible examples; see
#'   `data-raw/make_nc_tracts.R`.
#' @examples
#' data(nc_sids)
#' data(nc_sids_W)
#' summary(nc_sids$log_sids)
"nc_sids"

#' Queen-contiguity weights for the North Carolina SIDS counties
#'
#' A sparse, symmetric binary (`style = "B"`) queen-contiguity weights matrix over the
#' 100 North Carolina counties of [nc_sids], with row/column names equal to
#' `nc_sids$county`. Pass to [spmixqr()] via `spatial_W = spq_weights(nc_sids_W,
#' type = "supplied")` (or directly as `spatial_W`).
#'
#' @format A 100 x 100 sparse `dgCMatrix` (490 nonzero links, one connected component,
#'   no islands).
#' @source Queen contiguity (\code{spdep::poly2nb(queen = TRUE)}) of the \pkg{sf}
#'   North Carolina shapefile; built offline in `data-raw/make_nc_tracts.R`.
#' @examples
#' data(nc_sids_W)
#' dim(nc_sids_W)
"nc_sids_W"

#' Columbus, Ohio neighbourhood crime
#'
#' Residential crime in 49 neighbourhoods of Columbus, Ohio (1980), the canonical
#' spatial-econometrics dataset (Anselin 1988). Shipped as a tidy data frame with
#' neighbourhood centroids, to illustrate spatial-error (CAR) quantile regression on
#' contiguous areal units. The matching queen-contiguity weights are in [columbus_W].
#'
#' @format A data frame with 49 rows (neighbourhoods) and 6 variables:
#' \describe{
#'   \item{id}{Neighbourhood id (the spatial unit; matches `columbus_W` row/col order).}
#'   \item{crime}{Residential burglaries and vehicle thefts per 1000 households.}
#'   \item{income}{Household income (USD 1000).}
#'   \item{hoval}{Housing value (USD 1000).}
#'   \item{x, y}{Neighbourhood-polygon centroids (planar map units).}
#' }
#' @source The `columbus` data distributed with the \pkg{spData} package (Anselin,
#'   L. 1988, \emph{Spatial Econometrics: Methods and Models}). Built offline from the
#'   \pkg{spData} `columbus.gpkg` shapefile and redistributed as a built-in for
#'   reproducible examples; see `data-raw/make_columbus.R`.
#' @examples
#' data(columbus_crime)
#' summary(columbus_crime$crime)
"columbus_crime"

#' Queen-contiguity weights for the Columbus neighbourhoods
#'
#' A sparse, symmetric binary (`style = "B"`) queen-contiguity weights matrix over the
#' 49 Columbus neighbourhoods of [columbus_crime] (118 links, one connected component).
#' Pass to [spmixqr()] via `spatial_W = spq_weights(columbus_W, type = "supplied")`.
#'
#' @format A 49 x 49 sparse `dgCMatrix` (118 nonzero links).
#' @source Queen contiguity (\code{spdep::poly2nb(queen = TRUE)}) of the \pkg{spData}
#'   Columbus shapefile; built offline in `data-raw/make_columbus.R`.
#' @examples
#' data(columbus_W)
#' dim(columbus_W)
"columbus_W"

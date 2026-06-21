## Build the shipped Columbus crime dataset + queen-contiguity weights.
## Decoupled from spData/sf at package-build time (avoids the local sp segfault):
## run this once to produce data/columbus_crime.rda and data/columbus_W.rda.
suppressMessages({library(spData); library(spdep); library(sf); library(Matrix)})
shp <- system.file("shapes/columbus.gpkg", package = "spData")
col <- sf::st_read(shp, quiet = TRUE)
columbus_crime <- data.frame(
  id    = seq_len(nrow(col)),
  crime = col$CRIME,   # residential burglaries + vehicle thefts per 1000 households
  income = col$INC,    # household income ($1000)
  hoval = col$HOVAL,   # housing value ($1000)
  x = col$X, y = col$Y # polygon centroids (arbitrary planar units)
)
nb <- spdep::poly2nb(col, queen = TRUE)
W  <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)   # symmetric binary queen contiguity
columbus_W <- Matrix::Matrix(W, sparse = TRUE)
dimnames(columbus_W) <- list(NULL, NULL)
stopifmnot <- function(x) if (!x) stop("check failed")
stopifmnot(isSymmetric(unname(as.matrix(columbus_W))))
save(columbus_crime, file = "data/columbus_crime.rda", compress = "xz")
save(columbus_W,     file = "data/columbus_W.rda",     compress = "xz")
cat("saved columbus_crime (", nrow(columbus_crime), "x", ncol(columbus_crime),
    ") and columbus_W (", nrow(columbus_W), "links", sum(columbus_W)/2, ")\n")

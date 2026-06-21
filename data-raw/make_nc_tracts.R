## Build the shipped areal example from sf's North Carolina SIDS shapefile.
## Decoupled: build the data frame + queen-contiguity weights here, save .rda.
## Run in a clean R session (Rscript), never library(sp) interactively (segfault note).
suppressMessages({library(sf); library(spdep); library(Matrix); library(methods)})
nc <- st_read(system.file("shape/nc.shp", package="sf"), quiet=TRUE)
## queen contiguity nb on the 100 counties
nb <- poly2nb(nc, queen=TRUE)
W <- nb2mat(nb, style="B", zero.policy=TRUE)
ids <- as.character(nc$NAME)
dimnames(W) <- list(ids, ids)
## tidy data frame: SIDS rate, births, non-white births, region label
nc_sids <- data.frame(
  county   = ids,
  births   = as.numeric(nc$BIR74),
  sid      = as.numeric(nc$SID74),
  nwbirths = as.numeric(nc$NWBIR74),
  east     = as.numeric(st_coordinates(st_centroid(st_geometry(nc)))[,1]),
  north    = as.numeric(st_coordinates(st_centroid(st_geometry(nc)))[,2]),
  stringsAsFactors = FALSE
)
## response: log SIDS rate per 1000 births (add 0.5 continuity), predictor: prop non-white
nc_sids$sids_rate <- 1000 * (nc_sids$sid + 0.5) / (nc_sids$births + 1)
nc_sids$log_sids  <- log(nc_sids$sids_rate)
nc_sids$pnw       <- nc_sids$nwbirths / (nc_sids$births + 1)
## ship the weights as a plain sparse matrix (no spdep object needed at load)
nc_sids_W <- methods::as(methods::as(Matrix::Matrix(W, sparse=TRUE), "CsparseMatrix"),
                         "generalMatrix")
## dump to CSV for provenance, then save .rda
write.csv(nc_sids, "data-raw/nc_sids.csv", row.names=FALSE)
save(nc_sids, file="data/nc_sids.rda", compress="xz")
save(nc_sids_W, file="data/nc_sids_W.rda", compress="xz")
cat("nc_sids:", nrow(nc_sids), "counties;  W nnz:", length(nc_sids_W@x),
    " components:", spdep::n.comp.nb(nb)$nc, "\n")
cat("DONE\n")

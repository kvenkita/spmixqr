## Test environments
* Local: Windows 11, R 4.4.3
* GitHub Actions: ubuntu-latest (release), windows-latest, macOS-latest (planned)

## Submission
Version 0.2.0 — adds an optional per-regime conditional-autoregressive (CAR/ICAR)
spatial-error term (weights from queen/rook contiguity, distance, k-NN, or
user-supplied), with permutation Moran's I diagnostics. New Imports: spdep, Matrix,
methods (all used). spatial_error = FALSE reproduces 0.1.0.

## R CMD check results
0 errors | 1 warning (expected, see below) | NOTEs as below.

* New submission.
* Imports 'mixqr', which is currently distributed via GitHub
  (github.com/kvenkita/mixqr) and not yet on CRAN. `spmixqr` will be submitted only
  after 'mixqr' is accepted; until then a "Strong dependencies not in mainstream
  repositories" NOTE/WARNING is expected. This is the family submission order
  (mixqr -> spmixqr), not a defect.

## Downstream dependencies
None.

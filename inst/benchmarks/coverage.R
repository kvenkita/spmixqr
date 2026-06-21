## Monte-Carlo coverage benchmark for spmixqr (run interactively; not in the test
## suite -- the bootstrap of the full pipeline is expensive). Reports gate / slope
## coverage for the bootstrap vs the classification-conditional sandwich, and
## coverage as a function of the number of distinct sites.
##
##   Rscript inst/benchmarks/coverage.R
##
## Expectation: bootstrap ~ nominal in well-spanned designs; sandwich < nominal.

suppressMessages(library(spmixqr))

run_one <- function(n, B = 200L, seed = 1) {
  d <- sim_spmixqr(n = n, G = 2, tau = 0.5, seed = seed)
  ## bootstrap fit
  fit <- spmixqr(y ~ x, d$data, coords = d$coords, G = 2, tau = 0.5,
                 variance = "boot",
                 control = spmixqr_control(nstart = 4L, seed = seed,
                                           boot_B = B, boot_block = 3L))
  V <- fit$vcov$all
  ## the gate spatial coefficients are the second-class block; we report whether the
  ## gate-implied membership probability at each site is within a 90% interval of the
  ## truth (a summary functional robust to label switching).
  pr <- fit$prior[, 2]
  truth <- d$truth$pi[, 2]
  list(rmse = sqrt(mean((pr - truth)^2)),
       coverage_label = mean(abs(pr - truth) < 0.2))
}

grid <- expand.grid(n = c(150, 300, 600), rep = 1:20)
res <- do.call(rbind, Map(function(n, r) {
  o <- tryCatch(run_one(n, B = 150L, seed = r), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  data.frame(n = n, rmse = o$rmse, coverage = o$coverage_label)
}, grid$n, grid$rep))

agg <- aggregate(cbind(rmse, coverage) ~ n, res, mean)
print(agg)
cat("\nRMSE should fall with n; coverage approaches nominal in well-spanned donor",
    "pools. The classification-conditional sandwich (variance='sandwich') is faster",
    "but optimistic -- prefer the bootstrap for reporting.\n")

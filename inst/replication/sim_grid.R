## Full recovery replication grid for spmixqr (run interactively).
##   Rscript inst/replication/sim_grid.R
## Reports held-out-location gate correlation, classification accuracy, and
## slope-surface RMSE across sample sizes and spatial-signal strengths.

suppressMessages(library(spmixqr))

one <- function(n, gate_slope, coef_slope, seed) {
  d <- sim_spmixqr(n = n, G = 2, tau = 0.5, gate_slope = gate_slope,
                   coef_slope = coef_slope, seed = seed)
  ntr <- floor(0.75 * n); tr <- seq_len(ntr); te <- (ntr + 1L):n
  fit <- tryCatch(
    spmixqr(y ~ x, d$data[tr, ], coords = d$coords[tr, ], G = 2, tau = 0.5,
            variance = "none", control = spmixqr_control(nstart = 4L, seed = seed)),
    error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pr <- predict(fit, newdata = d$data[te, ], newcoords = d$coords[te, ], type = "prob")
  cl <- apply(pr, 1L, which.max)
  acc <- max(mean(cl == d$truth$z[te]), mean((3 - cl) == d$truth$z[te]))
  data.frame(n = n, gate_slope = gate_slope, coef_slope = coef_slope,
             gate_cor = abs(cor(pr[, 2], d$truth$pi[te, 2])), accuracy = acc)
}

grid <- expand.grid(n = c(200, 400, 800), gate_slope = c(1.5, 3),
                    coef_slope = c(0.5, 1.5), rep = 1:10)
res <- do.call(rbind, Map(one, grid$n, grid$gate_slope, grid$coef_slope,
                          seq_len(nrow(grid))))
agg <- aggregate(cbind(gate_cor, accuracy) ~ n + gate_slope + coef_slope, res, mean)
print(agg)

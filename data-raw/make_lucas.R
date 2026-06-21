## Build a shipped subsample of the Lucas County (Ohio) house-sales data for the
## mixture example. Decoupled from spData/sp (avoids the local sp segfault): run once.
suppressMessages({library(spData); library(sp)})
data(house, package = "spData")
df <- as.data.frame(house); cc <- sp::coordinates(house)
set.seed(7); ix <- sample(nrow(df), 2000)
lucas_house <- data.frame(
  price = df$price[ix],            # sale price (USD)
  tla   = df$TLA[ix],              # total living area (sq ft)
  age   = df$age[ix],              # age at sale (years)
  rooms = df$rooms[ix],
  x = cc[ix, 1], y = cc[ix, 2]     # projected coordinates (feet)
)
cat("age range:", range(lucas_house$age), " price:", range(lucas_house$price),
    " tla:", range(lucas_house$tla), "\n")
save(lucas_house, file = "data/lucas_house.rda", compress = "xz")
cat("saved lucas_house", nrow(lucas_house), "x", ncol(lucas_house), "\n")

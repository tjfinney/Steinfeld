# dimensions.R
# Get matrix dimensions.
# TJF 27/5/17

path <- "../dist"
regex <- "Rom-Steinfeld.+counts"
dist.mx.list <- list.files(path, pattern = regex)
lapply(dist.mx.list, function(x) {
  # Get dist object
  dist.mx.path <- sprintf("%s/%s", path, x)
  dist.mx <- read.csv(dist.mx.path, row.names=1)
  message(dist.mx.path, ": ", dim(dist.mx)[1], " x ", dim(dist.mx)[2])
})
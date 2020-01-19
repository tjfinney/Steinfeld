# DC-cut.R
# T. J. Finney, 2017-06-04
message("Cut DC tree to produce k groups...")

# Reset.
rm(list=ls(all=TRUE))

# Set parameters.
setwd('/home/tjf2n/Keep/www/tfinney.net/Views/restricted/Steinfeld/scripts')
message("Distance matrix (as CSV file): ", dist.path <- "../dist/Rom-Steinfeld-versions.15.SMD.csv")
message("Cutting height: ", height <- 0.5)

# Script

# Read distance matrix
dist.frame <- read.csv(dist.path, row.names=1)
colnames(dist.frame) <- rownames(dist.frame)

# Convert distance matrix to a distance object.
dist <- as.dist(dist.frame)

# Perform DC analysis.
require(cluster)
DC <- diana(dist)

# Cut tree at specified height
g <- cutree(as.hclust(DC), h=height)
# Get corresponding number of groups
k <- max(g)

# Print groups
cat("Group | Members\n")
cat("------|------------------------------------------------------------|------------\n")
lapply(1:k, function(n, groups) {
  cat(n, " | ", paste(names(groups)[groups==n], collapse=", "))
  cat("\n")
}, groups=g)

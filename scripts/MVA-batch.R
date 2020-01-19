# MVA-batch.R
# T. J. Finney, 2017-08-14
# Happy birthday Isa!
message("
This does the following kinds of analysis for distance matrices associated with a data matrix:
  * CMDS, DC, PAM for general distance matrices
  * ranked distances and group distances for special distance matrices.
Make distance matrices with \'dist-batch.R\' before running this script.
Run this from the R command line to avoid margin errors.
")

# Initialize
require(ape)
require(cluster)
require(rgl)
rm(list=ls(all=TRUE))
setwd('/home/tjf2n/Keep/www/tfinney.net/Views/restricted/Steinfeld/scripts')
message("Data matrix name: ", data.mx.name <- "Combined-Origen.csv")
message("No of groups for quantitative analysis: ", k <- 8)

# Functions

# CMDS
cmds.fn <- function(dist.mx.path, output.path) {
  dist.mx <- as.dist(read.csv(dist.mx.path, row.names=1))
  MDS <- cmdscale(dist.mx, k=3, eig=TRUE)
  # Calculate R-squared.
  MDS.dist <- dist(MDS$points, diag=TRUE, upper=TRUE)
  MDS.summary <- summary.lm(lm(dist.mx ~ MDS.dist))
  rsq <- MDS.summary$r.squared
  # Make result
  graphics.off()
  par3d(windowRect = c(0, 0, 600, 600))
  par3d(family = "bitmap")
  x <- MDS$points[,1]
  y <- MDS$points[,2]
  z <- MDS$points[,3]
  plot3d(x, y, z, xlab="axis 1", ylab="axis 2", zlab="axis 3", type='n', axes=TRUE, box=TRUE, sub=sprintf("R-squared = %0.2f", rsq))
  text3d(x, y, z, rownames(MDS$points), col=4)
  movie3d(spin3d(rpm=10), duration=6, fps=16, dir=".", movie=output.path)
}

# DC
dc.fn <- function(dist.mx.path, output.path) {
  dist.mx <- as.dist(read.csv(dist.mx.path, row.names=1))
  DC <- diana(dist.mx)
  graphics.off()
  par(bg="white")
  plot(DC, which=2, main="", xlab="", cex=0.8)
  dev.print(png, file=output.path, width=900, height=600)
}

# PAM
pam.fn <- function(dist.mx.path, output.path) {
  dist.mx <- as.dist(read.csv(dist.mx.path, row.names=1))
  # MSW plot
  N <- attr(dist.mx, "Size")
  x <- 1:(N - 1)
  y <- vector()
  for (n in x) { y[n] <- pam(dist.mx, n)$silinfo$avg.width }
  graphics.off()
  par(bg="white")
  plot(x, y, main=NULL, xlab="no. of groups", ylab="mean silhouette width")
  dev.print(png, file=output.path, width=600, height=600)
  # Table
  output.path <- gsub("png$", "txt", output.path)
  zz <- file(output.path, "w")
  cat("k | Groups | Poorly classified\n", file = zz)
  cat("------|------------------------------------------------------------|------------\n", file = zz)
  for (k in c(2,3,4,5,8,15,30)) {
    PAM <- pam(dist.mx, k)
    groups <- vector(mode = "character", length = k)
    for (n in 1:k) {
      group <- names(PAM$clustering[PAM$clustering == n])
      medoid <- PAM$medoids[n]
      group[group == medoid] <- sprintf("(%s)", medoid)
      groups[[n]] <- sprintf("[%s]", paste(group, collapse=" "))
    }
    groups <- paste(groups, collapse=" ")
    sil.widths <- sort(PAM$silinfo$widths[,3], decreasing=TRUE)
    poor <- paste(names(sil.widths[sil.widths < 0]), collapse=" ")
    cat(sprintf("%s | %s | %s\n", k, groups, poor), file = zz)
  }
  close(zz)
}

# Ranked distances
RD.fn <- function(refs, dist.mx.refs, counts.mx.refs, output.path) {
  zz <- file(output.path, "w")
  cat("Reference | Ranked distances\n", file = zz)
  cat("------|------------------------------------------------------------\n", file = zz)
  for (i in 1:length(refs)) {
    ref <- refs[i]
    dist.mx <- read.csv(dist.mx.refs[i], row.names=1)
    counts.mx <- read.csv(counts.mx.refs[i], row.names=1)
    p <- mean(as.dist(dist.mx))
    counts <- unlist(counts.mx[rownames(counts.mx)==ref])
    ref.fr <- dist.mx[rownames(dist.mx)==ref]
    ref.fr[2] <- round(ref.fr[1] * counts)
    ref.fr[3] <- counts
    ref.fr[4] <- qbinom(0.025, counts, p)/counts
    ref.fr[5] <- qbinom(0.975, counts, p)/counts
    ref.fr[6] <- (ref.fr[1] < ref.fr[4]) | (ref.fr[1] > ref.fr[5])
    ref.fr[7] <- row.names(ref.fr)
    ranked <- ref.fr[order(ref.fr[1]),]
    ranked <- subset(ranked, rownames(ranked) != ref)
    parts <- apply(ranked, 1, function(r) {
      # use gsub to drop spurious spaces
      test <- as.logical(gsub("[[:space:]]", "", r[6]))
      ratio <- gsub("[[:space:]]", "", sprintf("%s/%s", r[2], r[3]))
      sprintf("%s (%s)%s", r[7], ratio, if (test) "\\*" else "")
    })
    cat(sprintf("%s | %s\n", ref, paste(parts, collapse="; ")), file = zz)
  }
  close(zz)
}

# Grouped distances
GD.fn <- function(dist.mx.def, k, refs, dist.mx.refs, output.path) {
  def.mx <- as.dist(read.csv(dist.mx.def, row.names=1))
  # Do PAM for with k groups
  def.pam <- pam(def.mx, k)
  medoids <- def.pam$medoids
  gr.fn <- function(m, p) {
    i <- which(p$medoids == m)
    names(p$clustering[p$clustering == i])
  }
  groups <- sapply(medoids, gr.fn, def.pam)
  row.fn <- function(n, rr, mm, gg) {
    ref <- rr[n]
    ref.mx <- read.csv(mm[n], row.names=1)
    dd <- ref.mx[rownames(ref.mx)==ref]
    sapply(gg, function(g, distances) {
      vals <- subset(distances, rownames(distances) %in% g)
      avg <- sprintf("%.3f", round(mean(unlist(vals)), 3))
      avg <- gsub("NaN", "NA", avg)
    }, dd)
  }
  rows <- sapply(1:length(refs), row.fn, refs, dist.mx.refs, groups)
  # Write output
  zz <- file(output.path, "w")
  cat(paste(c("Reference", medoids), collapse = " | "), "\n", file = zz)
  cat(paste(rep("--------", k + 1), collapse = "|"), "\n", file = zz)
  sapply(1:length(refs), function(i) {
    parts <- c(refs[i], rows[,i])
    cat(paste(parts, collapse = " | "), "\n", file = zz)
  })
  close(zz)
}

# Script

# Check for required directories
if (!(
  file_test("-d", "../dist") &&
  file_test("-d", "../cmds") &&
  file_test("-d", "../dc") &&
  file_test("-d", "../pam") &&
  file_test("-d", "../RD") &&
  file_test("-d", "../GD")
))
stop("Directories dist, cmds, dc, pam, RD, and GD must exist as siblings of the directory containing this script.", call.=FALSE)

# Make matrix paths
# General distance matrix
data.mx.prefix <- sub(".csv$", "", data.mx.name)
regex <- sprintf("%s\\.\\d+\\.SMD", data.mx.prefix)
dist.mx.def <- list.files("../dist", pattern = regex)
# Special (i.e. witness-specific) distance matrices, counts matrices, and names
regex <- sprintf("%s.+SMD", data.mx.prefix)
dist.mx.refs <- list.files("../dist", pattern = regex)
dist.mx.refs <- dist.mx.refs[dist.mx.refs != dist.mx.def]
counts.mx.refs <- gsub("SMD", "counts", dist.mx.refs)
refs <- sub("^[^.]+.([^.]+).+", "\\1", counts.mx.refs)
# Add prefix
dist.mx.def <- sprintf("../dist/%s", dist.mx.def)
dist.mx.refs <- sprintf("../dist/%s", dist.mx.refs)
counts.mx.refs <- sprintf("../dist/%s", counts.mx.refs)

# Do MVA (CMDS, DC, PAM) for general distance matrix
output.path = gsub(".csv$", "", gsub("/dist/", "/cmds/", dist.mx.def))
cmds.fn(dist.mx.def, output.path)
output.path = gsub(".csv$", ".png", gsub("/dist/", "/dc/", dist.mx.def))
dc.fn(dist.mx.def, output.path)
output.path = gsub(".csv$", ".png", gsub("/dist/", "/pam/", dist.mx.def))
pam.fn(dist.mx.def, output.path)

# Do ranked distances and group distances for special witnesses
output.path = gsub("\\.\\d+.SMD.csv$", ".txt", gsub("/dist/", "/RD/", dist.mx.def))
RD.fn(refs, dist.mx.refs, counts.mx.refs, output.path)
output.path = gsub("/RD/", "/GD/", output.path)
GD.fn(dist.mx.def, k, refs, dist.mx.refs, output.path)

# dist-batch.R
# T. J. Finney, 2017-05-14
# This incorporates code suggested by Bill Venables.

message("
Make distance and count matrices from a data matrix.
One general and a number of special distance matrices are made.
The minimum number of sampling points for the general one is 15;
the min. no. for special distance matrices is set by script variable \'samp\'.
A special distance matrix is made for each witness that is not in the general distance matrix but has more than \'samp\' sampling points.
")

# Initialize
require(cluster)
rm(list=ls(all=TRUE))
setwd('/home/tjf2n/Keep/www/tfinney.net/Views/restricted/Steinfeld/scripts')
message("Data matrix name: ", data.mx.name <- "Combined-Origen.csv")
message("Min. no. of sampling points for special distance matrices: ", samp <- 8)

# Functions

# Make path name
# Arguments:
# * data.mx.path: path to data matrix
# * ref: item to keep (if possible); use "" if no ref. item
# * limit: minimum number of pairwise defined states
# * type: type descriptor (e.g. 'SMD', 'counts')
# Return:
# * path name as a string.

path.name.fn <- function(data.mx.path, ref, limit, type) {
  sub(
    ".csv$",
    if (ref == "") sprintf(".%s.%s.csv", limit, type)
    else sprintf(".%s.%s.%s.csv", ref, limit, type),
    sub("/data/", "/dist/", data.mx.path)
  )
}

# Make distance and counts matrices from a data matrix:
# * items are eliminated until all numbers of pairwise defined states exceed a specified limit
# Arguments:
# * data.mx.path: path to data matrix (a CSV file comprised of nominal data)
# * ref: item to keep (if possible); use "" if no ref. item
# * limit: minimum number of pairwise defined states
# * quiet: whether to be quiet (i.e. not generate messages)
# Return:
# * items in the distance and counts matrices as a vector.
# Note:
# The distance and count matrices are made as side-effects.

dist_count.fn <- function(data.mx.path, ref, limit, quiet) {
  # Read input
  input.fr <- read.csv(data.mx.path, row.names=1, colClasses="factor")
  # Make pairwise defined counts matrix
  input.mx <- data.matrix(input.fr)
  input.mx[] <- !is.na(input.mx)
  counts.mx <- tcrossprod(input.mx)
  # Stage 1: drop objects without enough defined states
  drop <- (diag(counts.mx) < limit)
  if (!quiet) {
    message(ref)
    message("drop items with < ", limit, " defined states:")
    message(paste(rownames(input.fr)[drop], collapse=" "))
  }
  st1.fr <- input.fr[!drop,]
  counts.mx <- counts.mx[!drop, !drop]
  # Stage 2: drop items without enough pairwise defined states
  while((m <- min(counts.mx)) < limit) {
    # List items with least counts
    least <- which(counts.mx == m, arr.ind=TRUE)[,1]
    # Choose worst defined item which is not ref
    wd <- sort(diag(counts.mx)[least[names(least) != ref]])[1]
    i <- which(names(least) == names(wd))[1]
    # Drop it
    counts.mx <- counts.mx[-least[i], -least[i]]
  }
  drop <- !(rownames(st1.fr) %in% rownames(counts.mx))
  if (!quiet) {
    message("drop items causing < ", limit, " pairwise defined states:")
    message(paste(rownames(st1.fr)[drop], collapse=" "))
  }
  st2.fr <- st1.fr[!drop,]
  # Min. count
  if (!quiet) {
    message("Min. no. of sampling points: ", min(counts.mx))
  }
  # Make outputs
  output.dist.mx <- as.matrix(daisy(st2.fr))
  # Write matrices
  write.csv(round(output.dist.mx, digits=3), path.name.fn(data.mx.path, ref, limit, "SMD"))
  write.csv(counts.mx, path.name.fn(data.mx.path, ref, limit, "counts"))
  # Return list of items in matrices
  rownames(counts.mx)
}

# Script

# Check for required directories
if (!(
  file_test("-d", "../data") &&
  file_test("-d", "../dist")
))
stop("Directories data and dist must exist as siblings of the directory containing this script.", call.=FALSE)

# Make path for data matrix
data.mx.path <- sprintf("../data/%s", data.mx.name)

# Check data matrix path
if (!file_test("-f", data.mx.path))
stop("Data matrix ", data.mx.path, " does not exist.", call.=FALSE)
if (length(grep(".csv$", data.mx.path)) == 0)
stop("Data matrix name must have \'.csv\' extension", call.=FALSE)

# Make general distance and count matrices

message("Make general distance and count matrices...")
general.names <- dist_count.fn(data.mx.path, "", 15, FALSE)

# Make special distance and count matrices
input.fr <- read.csv(data.mx.path, row.names=1, colClasses="factor")
input.mx <- data.matrix(input.fr)
input.mx[] <- !is.na(input.mx)
counts.mx <- tcrossprod(input.mx)
enough <- names(diag(counts.mx)[diag(counts.mx) >= samp])
good <- setdiff(enough, general.names)
message("Make distance and count matrices for ", paste(good, collapse=" "))
for (ref in good) {
  dist_count.fn(data.mx.path, ref, samp, TRUE)
}

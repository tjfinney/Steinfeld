setwd('/home/tjf2n/Keep/www/tfinney.net/Views/restricted/Steinfeld/scripts')
A <- "../dist/Cor-B-Steinfeld-Origen.15.SMD.csv"
B <- "../dist/Cor-B-Mallett-versions.15.SMD.csv"
A.fr <- read.csv(A, row.names=1, colClasses="factor")
B.fr <- read.csv(B, row.names=1, colClasses="factor")
message(sprintf("In %s and %s:", A, B))
message(cat(intersect(rownames(A.fr), rownames(B.fr))))
message(sprintf("In %s but not %s:", A, B))
message(cat(setdiff(rownames(A.fr), rownames(B.fr))))
message(sprintf("In %s but not %s:", f2, f1))
message(cat(setdiff(rownames(B.fr), rownames(A.fr))))

# Draft analysis of motifs in topics k = 1 through 10.
homer <- readRDS(file.path("../output/Buenrostro_2018",
                           "binarized/filtered_peaks",
                           "homer-buenrostro2018-k=10-noshrink.rds"))

# Compile the motif enrichment p-values into a single table, after
# removing duplicate results.
motifs <- sort(unique(homer$k1[,"Motif Name"]))
n      <- length(motifs)
pvals  <- matrix(0,n,k)
rownames(pvals) <- motifs
colnames(pvals) <- paste0("k",1:k)
for (i in 1:k) {
  dat  <- homer[[i]]
  rows <- which(!duplicated(dat[,"Motif Name"]))
  dat  <- dat[rows,]
  rownames(dat) <- dat[,"Motif Name"]
  pvals[,i] <- dat[motifs,"P-value"]
}
pvals <- as.data.frame(pvals)


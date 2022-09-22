library(fastTopics)
load("../data/cicero_gene_kidney.RData")
load(file.path("../output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))
peaks  <- rownames(de$postmean)
counts <- rep(0,length(peaks))
names(counts) <- peaks
genes  <- names(cicero_gene)
n      <- length(genes)
for (i in 1:n) {
  cat(i,"")
  gene      <- genes[i]    
  i         <- intersect(cicero_gene[[gene]],peaks)
  counts[i] <- counts[i] + 1
}
cat("\n")

# In this short script we use the results of the gene-by-gene adaptive
# shrinkage in compute_gene_scores_cusanovich2018_kidney_k10.R to
# generate a de_analysis data structure for genes instead of peaks.
library(fastTopics)
library(tools)
load(file.path("../output/Cusanovich_2018/tissues",
               "gene-scores-cusanovich2018-kidney-k=10.RData"))
k <- 10
n <- length(gene_scores)
A <- matrix(as.numeric(NA),n,k)
genes       <- names(gene_scores)
rownames(A) <- genes
colnames(A) <- paste0("k",1:k)
de_gene <- list(postmean = A,se = A,z = A,lfsr = A)
class(de_gene) <- c("topic_model_de_analysis","list")
for (i in 1:n) {
  cat(sprintf("%d (%s) ",i,genes[i]))
  res <- gene_scores[[i]]
  if (is.list(res)) {
    for (j in 1:k) {
      row <- which.min(res$lfsr[,j])
      de_gene$postmean[i,j] <- res$b[row,j]
      de_gene$se[i,j]       <- res$se[row,j]
      de_gene$z[i,j]        <- res$z[row,j]
      de_gene$lfsr[i,j]     <- res$lfsr[row,j]
    }
  }
}
cat("\n")
save(list = "de_gene",file = "de-gene-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("de-gene-cusanovich2018-kidney-k=10.RData")

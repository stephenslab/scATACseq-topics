# In this short script we use the results of the gene-by-gene adaptive
# shrinkage in compute_gene_scores_cusanovich2018_kidney_k10.R to
# generate a de_analysis data structure for genes instead of peaks.
library(fastTopics)
library(tools)
load("gene-scores-cusanovich2018-kidney-k=10.RData")
n <- length(gene_scores)
k <- length(gene_scores[[1]]$logLR)
A <- matrix(as.numeric(NA),n,k)
genes       <- names(gene_scores)
rownames(A) <- genes
colnames(A) <- paste0("k",1:k)
de_gene <- list(coef = A,logLR = A)
for (i in 1:n) {
  cat(i,"")
  res <- gene_scores[[i]]
  for (j in 1:k) {
    de_gene$coef[i,j]  <- res$coef[j]
    de_gene$logLR[i,j] <- res$logLR[j]
  }
}
cat("\n")
save(list = "de_gene",file = "de-gene-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("de-gene-cusanovich2018-kidney-k=10.RData")

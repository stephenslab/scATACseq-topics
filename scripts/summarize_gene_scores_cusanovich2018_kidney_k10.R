# In this short script we use the results of the gene-by-gene adaptive
# shrinkage in compute_gene_scores_cusanovich2018_kidney_k10.R to
# generate a de_analysis data structure for genes instead of peaks.
library(fastTopics)
library(tools)
load(file.path("../output/Cusanovich_2018/tissues",
               "gene-scores-cusanovich2018-kidney-k=10.RData"))
n <- length(gene_scores)
k <- length(gene_scores[[1]]$logLR)
A <- matrix(as.numeric(NA),n,k)
genes       <- names(gene_scores)
rownames(A) <- genes
colnames(A) <- paste0("k",1:k)
de_gene <- list(postmean = A,se = A,z = A,lfsr = A)
class(de_gene) <- c("topic_model_de_analysis","list")
for (i in 1:n) {
  cat(i,"")
  res <- gene_scores[[i]]
  for (j in 1:k) {
    r <- which.min(res$lfsr[,j])

    # To simplify creation of the volcano plots, here we store the
    # mean coefficients (coef) in the "postmean" slot, and we store
    # the log-likelihood ratios in the "z" slow. The standard errors
    # (se) and lfsr's are filled in as well, although they are less
    # essential.
    de_gene$postmean[i,j] <- res$coef[j]
    de_gene$z[i,j]        <- res$logLR[j]
    de_gene$se[i,j]       <- res$se[r,j]
    de_gene$lfsr[i,j]     <- res$lfsr[r,j]
  }
}
cat("\n")
save(list = "de_gene",file = "de-gene-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("de-gene-cusanovich2018-kidney-k=10.RData")

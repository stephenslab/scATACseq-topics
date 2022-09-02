# TO DO: Explain here what this script is for, and how to use it.
library(fastTopics)
library(tools)
load(file.path("../output/Cusanovich_2018/tissues",
               "gene-scores-cusanovich2018-kidney-k=10.RData"))
k <- 10
n <- length(gene_scores)
A <- matrix(as.numeric(NA),n,k)
rownames(A) <- names(gene_scores)
colnames(A) <- paste0("k",1:k)
de_gene <- list(postmean = A,se = A,z = A,lfsr = A)
class(de) <- c("topic_model_de_analysis","list")

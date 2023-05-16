# Here we perform adaptive shrinkage separately for the peaks near
# each gene in order to produce gene-wise statistics from the results
# at the level of chromatin accessibility peaks. For this analysis we
# use the de_analysis results generated from the
# de_analysis_cusanovich2018_kidney_k10.R script.
library(fastTopics)
library(ashr)
library(tools)
source("../code/ash.R")

# The "prior sample size" used when running ash on the peaks linked to
# each gene.
n0 <- 20

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the sparse matrix associating genes with peaks.
load("../data/cicero_gene_kidney.RData")

# Load the results of the differential accessibility analysis using
# the k = 10 topic model.
load("de-cusanovich2018-kidney-k=10-noshrink.RData")

# Run ash once on all peaks and all topics.
b  <- as.vector(de$postmean)
z  <- as.vector(de$z)
se <- b/z
res0 <- ash(b,se,mixcompdist = "normal",method = "shrink",outputlevel = 1)
g0 <- res0$fitted_g

# For each gene, perform adaptive shrinkage on the de_analysis results
# for all peaks near the gene. The adaptive shrinkage is performed
# separately for each topic.
k     <- ncol(de$postmean)
peaks <- rownames(de$postmean)
genes <- colnames(cicero_gene_mat)
gene_scores <- vector("list",length(genes))
names(gene_scores) <- genes
for (gene in genes) {
  cat(gene,"")

  # Get the peaks near the gene.
  rows <- which(cicero_gene_mat[,gene] > 0)
  rows <- which(is.element(peaks,rownames(cicero_gene_mat)[rows]))
  n    <- length(rows)
  if (n > 0) {

    # Set up storage for the ash outputs.
    res <- list(logLR = rep(as.numeric(NA),k),
                coef  = rep(0,k),
                b     = matrix(0,n,k),
                se    = matrix(as.numeric(NA),n,k),
                z     = matrix(0,n,k),
                lfsr  = matrix(1,n,k))
    rownames(res$b)    <- rownames(de$postmean)[rows]
    rownames(res$se)   <- rownames(de$postmean)[rows]
    rownames(res$z)    <- rownames(de$postmean)[rows]
    rownames(res$lfsr) <- rownames(de$postmean)[rows]
    colnames(res$b)    <- paste0("k",1:k)
    colnames(res$se)   <- paste0("k",1:k)
    colnames(res$z)    <- paste0("k",1:k)
    colnames(res$lfsr) <- paste0("k",1:k)
    names(res$coef)  <- paste0("k",1:k)
    names(res$logLR) <- paste0("k",1:k)
    
    # Perform the adaptive shrinkage separately for each topic.
    for (i in 1:k) {

      # Set up the ash inputs.
      b  <- de$postmean[rows,i]
      z  <- de$z[rows,i]
      se <- b/z
      se[z == 0] <- as.numeric(NA)
      se[b == 0] <- 0

      # Perform the adaptive shrinkage step.
      if (any(!is.na(se))) {
        out <- ash_test_enrich(b,se,g0,prior = 1.01 + n0*g0$pi)

        # Store the ash results.
        res$b[,i]    <- out$b
        res$se[,i]   <- out$se
        res$z[,i]    <- out$z
        res$lfsr[,i] <- out$lfsr
        res$coef[i]  <- out$coef
        res$logLR[i] <- out$logLR
      }
    }
    
    # Store the results for that gene.
    gene_scores[[gene]] <- res
  }
}
cat("\n")

# Remove genes that do not have any peaks with differentially
# accessibility results.
i <- which(!sapply(gene_scores,is.null))
gene_scores <- gene_scores[i]

# Save the results to an .RData file.
save(list = "gene_scores",
     file = "gene-scores-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("gene-scores-cusanovich2018-kidney-k=10.RData")

# Here we perform adaptive shrinkage separately for the peaks near
# each gene in order to produce gene-wise statistics from the results
# at the level of chromatin accessibility peaks. For this analysis we
# use the de_analysis results generated from the
# de_analysis_cusanovich2018_kidney_k10.R script.
library(fastTopics)
library(ashr)
library(tools)
source("../code/ash.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the sparse matrix associating genes with peaks.
load("../data/cicero_gene_kidney.RData")

# Load the results of the differential accessibility analysis using
# the k = 10 topic model.
load(file.path("../output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))

# Run ash once on all peaks and all topics.
b  <- as.vector(de$postmean)
z  <- as.vector(de$z)
se <- b/z
res0 <- ash(b,se,mixcompdist = "normal",method = "shrink",outputlevel = 1)

# For each gene, perform adaptive shrinkage on the de_analysis results
# for all peaks near the gene.
peaks <- rownames(de$postmean)
genes <- colnames(cicero_gene_mat)
gene_scores <- vector("list",length(genes))
names(gene_scores) <- genes
for (gene in genes) {
  cat(gene,"")

  # Get the peaks near the gene.
  rows <- which(cicero_gene_mat[,gene] > 0)
  rows <- which(is.element(peaks,rownames(cicero_gene_mat)[rows]))
  if (length(rows) > 0) {

    # Set up the ash inputs.
    b  <- de$postmean[rows,,drop = FALSE]
    z  <- de$z[rows,,drop = FALSE]
    se <- b/z
    se[z == 0] <- as.numeric(NA)
    se[b == 0] <- 0

    # Perform adaptive shrinkage.
    res <- shrink_estimates(b,se,g = res0$fitted_g,fixg = length(rows) < 10)
    res$loglik <- res$ash$loglik
    res$logLR  <- res$ash$logLR
    
    # Store the ash results.
    gene_scores[[gene]] <- res[c("loglik","logLR","b","se","z","lfsr")]
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

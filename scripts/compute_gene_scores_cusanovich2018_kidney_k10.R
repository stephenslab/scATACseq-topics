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

# Load the data structure associating genes with peaks.
load("../data/cicero_gene.RData")

# Load the results of the differential accessibility analysis using
# the k = 10 topic model.
load(file.path("../output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))

# Run ash once on all peaks for a selected topic to obtain a default
# prior. Here we use topic k = 8
# because it seems to correspond well to the Loop of Henle (LoH) cell
# type.
j <- "k8"
b <- de$postmean[,j]
z <- de$z[,j]
se <- b/z
res0 <- ash(b,se,mixcompdist = "normal",method = "shrink")

# res1 <- ash(b[1:100],se[1:100],mixcompdist = "normal",method = "shrink",
#             g = res$fitted_g,fixg = TRUE)
# res2 <- ash(b[1:100],se[1:100],mixcompdist = "normal",method = "shrink",
#             g = res$fitted_g)

# For each gene, perform adaptive shrinkage on the de_analysis
# results for all peaks near the gene.
peaks <- rownames(de$postmean)
genes <- names(cicero_gene)
gene_scores <- vector("list",length(genes))
names(gene_scores) <- genes
for (gene in genes) {
  cat("%s ",gene)
  
  # Get the peaks near the gene.
  rows <- which(is.element(peaks,cicero_gene[[gene]]))
  
  # Set up the ash inputs.
  b  <- de$postmean[rows,]
  z  <- de$z[rows,]
  se <- b/z
  se[z == 0] <- as.numeric(NA)
  se[b == 0] <- 0

  # Perform adaptive shrinkage.
  res <- shrink_estimates(b,se,g = res0$fitted_g,fixg = length(rows) < 5)

  # Store the ash results.
  gene_scores[[gene]] <- res[c("b","se","z","lfsr")]
}
cat("\n")

# Save the results to an .RData file.
save(list = "gene_scores",
     file = "gene-scores-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("gene-scores-cusanovich2018-kidney-k=10.RData")

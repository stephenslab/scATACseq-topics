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

# Load the base-pair positions of the genes for the mm9 Mouse Genome
# Assembly.
load("../data/mm9_seq_gene.RData")
chromosomes <- c(1:19,"X","Y")
seq_gene <- transform(seq_gene,chromosome = as.character(chromosome))
seq_gene <- subset(seq_gene,is.element(chromosome,chromosomes))
seq_gene <- transform(seq_gene,chromosome = factor(chromosome,chromosomes))

# Load the results of the differential accessibility analysis using
# the k = 10 topic model.
load(file.path("../output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))

# Get the base-pair positions of the peaks.
feature_names <- rownames(de$postmean)
out           <- strsplit(feature_names,"_")
positions     <- data.frame(chr   = sapply(out,"[[",1),
                            start = sapply(out,"[[",2),
                            end   = sapply(out,"[[",3),
                            name  = feature_names,
                            stringsAsFactors = FALSE)
positions <- transform(positions,
                       chr   = factor(substr(chr,4,100),chromosomes),
                       start = as.numeric(start),
                       end   = as.numeric(end))

# Run ash once on all peaks for a selected topic to obtain a
# reasonable setting for "mixsd". Here we use topic k = 8 because it
# seems to correspond well to the Loop of Henle (LoH) cell type.
j <- "k8"
b <- de$postmean[,j]
z <- de$z[,j]
se <- b/z
res <- ash(b,se,mixcompdist = "normal",method = "shrink")
mixsd <- res$fitted_g$sd

# For each gene, and for each topic, perform adaptive shrinkage (ash)
# on the de_analysis results for all peaks near the gene.
d <- 2e5
n <- nrow(seq_gene)
k <- ncol(de$postmean)
gene_scores <- rep(list(NA),n)
names(gene_scores) <- seq_gene$feature_name
for (i in 1:n) {

  # Get the peaks near the gene.
  gene_dat <- seq_gene[i,]
  cat(sprintf("%d (%s) ",i,gene_dat$feature_name))
  rows <- which(with(positions,
                     chr   == gene_dat$chromosome &
                     start >= gene_dat$chr_start - d &
                     end   <= gene_dat$chr_stop + d))
  if (length(rows) > 1) {

    # Set up the ash inputs.
    b  <- de$postmean[rows,]
    z  <- de$z[rows,]
    se <- b/z
    se[z == 0] <- as.numeric(NA)
    se[b == 0] <- 0

    # Perform adaptive shrinkage separately for each topic.
    res <- shrink_estimates(b,se,mixsd)

    # Store the ash results.
    rownames(res$b)    <- positions[rows,"name"]
    rownames(res$se)   <- positions[rows,"name"]
    rownames(res$lfsr) <- positions[rows,"name"]
    colnames(res$b)    <- paste0("k",1:k)
    colnames(res$se)   <- paste0("k",1:k)
    colnames(res$lfsr) <- paste0("k",1:k)
    gene_scores[[i]]   <- res[c("b","se","z","lfsr")]
  }
}
cat("\n")

# Save the results to an .RData file.
save(list = "gene_scores",
     file = "gene-scores-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("gene-scores-cusanovich2018-kidney-k=10.RData")

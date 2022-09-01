# TO DO: Explain here what this script does, and how to use it.
library(fastTopics)
library(ashr)
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

# TO DO: Run ash once on all peaks for a selected topic to obtain a
# reasonable setting for mixsd.

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
  m <- length(rows)
  if (m > 1) {

    # If there are fewer than 20 peaks near the gene, do not estimate
    # the ash prior.
    fixg <- (m >= 20)

    # Set up the ash inputs.
    b  <- de$postmean[rows,]
    z  <- de$z[rows,]
    se <- b/z
    se[z == 0] <- as.numeric(NA)
    se[b == 0] <- 0

    # Perform adaptive shrinkage separately for each topic.
    #
    # TO DO: Add a prior to encourage non-sparse weights.
    #
    res <- fastTopics:::shrink_estimates(b,se,prior = rep(1,16))

    # Store the ash results.
    rownames(res$b)    <- positions[rows,"name"]
    rownames(res$se)   <- positions[rows,"name"]
    rownames(res$lfsr) <- positions[rows,"name"]
    colnames(res$b)    <- paste0("k",1:k)
    colnames(res$se)   <- paste0("k",1:k)
    colnames(res$lfsr) <- paste0("k",1:k)
    gene_scores[[i]] <- out[c("b","se","lfsr")]
  }
}
cat("\n")

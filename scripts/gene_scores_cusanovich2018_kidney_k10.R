# TO DO: Explain here what this script does, and how to use it.
library(fastTopics)
library(ashr)

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

# For each gene, and for each topic, perform adaptive shrinkage (ash)
# on the de_analysis results for all peaks near the gene.
d <- 2e5
n <- nrow(seq_gene)
k <- ncol(de$postmean)
gene_scores <- vector("list",n)
numpeaks <- rep(0,n)
names(gene_scores) <- seq_gene$feature_name
names(numpeaks) <- seq_gene$feature_name
for (i in 1:n) {
  gene_dat <- seq_gene[i,]
  cat(sprintf("%d (%s) ",i,gene_dat$feature_name))
  rows <- which(with(positions,
                     chr   == gene_dat$chromosome &
                     start >= gene_dat$chr_start - d &
                     end   <= gene_dat$chr_stop + d))
  m    <- length(rows)
  numpeaks[i] <- m
  if (m > 1) {
    res <- list(b      = matrix(0,m,k),
                se     = matrix(0,m,k),
                lfsr   = matrix(0,m,k),
                svalue = matrix(0,m,k))
    ash_dat <- data.frame(b =  de$postmean
  # for (j in 1:k) {
  se <- with(out,postmean/z)
    se[out$z == 0] <- as.numeric(NA)
    se[out$postmean == 0] <- 0
    res          <- shrink_estimates(out$postmean,se,...)
  #   dat <- data.frame(b = de$postmean[rows,k])
  #   dat <- transform(dat,se = postmean/z)
  #   fit <- ash(dat$b,pdat$se,mixcompdist = "normal",method = "shrink")
}
cat("\n")


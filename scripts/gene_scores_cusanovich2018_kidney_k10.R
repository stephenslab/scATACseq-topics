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

# For each gene, perform adaptive shrinkage (ash) on the de_analysis
# results for all peaks near the gene.
d <- 2e5
n <- nrow(seq_gene)
gene_scores <- vector("list",n)
names(gene_scores) <- seq_gene$feature_name
for (i in 1:n) {
  cat(sprintf("%d (%s) ",i,seq_gene[i,"feature_name"]))
}
cat("\n")

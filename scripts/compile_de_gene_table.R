# A short script to compile a CSV file containing the results of the
# gene enrichment analysis generated using scripts
# compute_gene_scores_cusanovich2018_kidney_k10.R and
# summarize_gene_scores_cusanovich2018_kidney_k10.R.
library(pathways)
data(gene_sets_mouse)

# Load the results of the gene enrichment analysis.
load(file.path("../output/Cusanovich_2018/tissues",
               "de-gene-cusanovich2018-kidney-k=10.RData"))

# Align the gene_info and de_gene data structures.
gene_info     <- gene_sets_mouse$gene_info
gene_info     <- subset(gene_info,!is.na(Ensembl))
gene_info     <- subset(gene_info,!duplicated(Ensembl))
ids           <- gene_info$Ensembl
ids           <- intersect(ids,rownames(de_gene$coef))
de_gene$coef  <- de_gene$coef[ids,]
de_gene$logLR <- de_gene$logLR[ids,]
rows          <- match(rownames(de_gene$coef),gene_info$Ensembl)
gene_info     <- gene_info[rows,]

# Compile the enrichment results for all topics into a single table.
k   <- ncol(de_gene$coef)
dat <- NULL
for (i in 1:k) {
  x <- data.frame(topic   = i,
                  ensembl = gene_info$Ensembl,
                  symbol  = gene_info$Symbol,
                  coef    = de_gene$coef[,i],
                  logLR   = de_gene$logLR[,i],
                  stringsAsFactors = FALSE)
  dat <- rbind(dat,x)
}

# Filter out genes with logLR <= 20.
dat <- subset(dat,logLR > 20)

# Reorder the genes by topic, then by logLR.
rows <- with(dat,order(topic,-logLR))
dat  <- dat[rows,]
rownames(dat) <- NULL

# Write the data frame to a CSV file.
dat <-
  transform(dat,
    topic = paste0("k",topic),
    coef  = format(round(coef,digits = 3),trim = TRUE,scientific = FALSE),
    logLR = format(round(logLR,digits = 3),trim = TRUE,scientific = FALSE))
write.csv(dat,"gene_enrich_cusanovich2018_kidney_k10.csv",
          quote = FALSE,row.names = FALSE)

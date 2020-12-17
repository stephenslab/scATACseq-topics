#! /usr/bin/env Rscript
# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# compute gene scores based on the weighted average of region scores,
# and perform gene-set enrichment analysis based on a multinomial topic model.

library(optparse)
library(tools)
library(Matrix)
library(fastTopics)
library(fgsea)
library(dplyr)
library(tidyr)
library(GenomicRanges)
source("../code/gsea.R")
source("../code/gene_annotation.R")
source("../code/gene_scores.R")

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--modelfit",type = "character",default="fit.rds")
parser <- add_option(parser,"--genesets",type="character",default="gene_sets.RData")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.RData")
parser <- add_option(parser,c("--genescore","-s"),type = "character",default = "genebody")
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=1000)
parser <- add_option(parser,"--extrapolate",action = "store_true")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countsfile  <- out$counts
prefitfile  <- out$prefit
outfile     <- out$out
method      <- out$method
numiter     <- out$numiter
extrapolate <- !is.null(out$extrapolate)
nc          <- out$nc
rm(parser,out)


#The ATAC-seq data was from the `mm9` version of mouse genome, so we load the TxDb and OrgDb for mouse `mm9` from Bioconductor.
suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm9.knownGene))
suppressPackageStartupMessages(library(org.Mm.eg.db))
TxDb  <- TxDb.Mmusculus.UCSC.mm9.knownGene
OrgDb <- org.Mm.eg.db

# SCRIPT SETTINGS
# ---------------
genesetfile  <- "/project2/mstephens/kevinluo/GSEA/pathways/output/gene_sets_mouse.RData"
countsfile   <- "/project2/mstephens/kevinluo/scATACseq-topics/data/Cusanovich_2018/processed_data/Cusanovich_2018.RData"
modelfitfile <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018/fit-Cusanovich2018-scd-ex-k=13.rds"
out.dir      <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018"
outfile      <- file.path(out.dir, "postfit-Cusanovich2018-scd-ex-k=13.RData")

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# LOAD MODEL FIT
# --------------
cat(sprintf("Loading Poisson NMF model fit from %s\n",modelfitfile))
fit <- readRDS(modelfitfile)$fit

# COMPUTE REGION Z-SCORES
# -----------------------
# Perform differential accessbility analysis using the multinomial topic model.
cat("Computing differential accessbility statistics from topic model.\n")
timing <- system.time(diff_count_res <- diff_count_analysis(fit,counts))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
saveRDS(diff_count_res, file.path(out.dir, "diff-count-Cusanovich2018-13topics.rds"))
# diff_count_res <- readRDS(file.path(out.dir, "diff-count-Cusanovich2018-13topics.rds"))

# COMPUTE GENE SCORES
# -------------------
# Extract genomic coordinates for ATAC-seq regions
region_scores <- diff_count_res$Z
regions <- data.frame(x = rownames(region_scores)) %>% separate(x, c("chr", "start", "end"), "_")
regions <- regions %>% mutate_at(c("start", "end"), as.numeric)

# Load gene annotations
genes <- GenomicFeatures::genes(TxDb)
genes <- as.data.frame(genes)
genes <- genes[,c("seqnames", "start", "end", "strand", "gene_id")]
colnames(genes) <- c("chr", "start", "end", "strand", "gene_id")

# Compute the gene scores
gene_scores <- compute_gene_scores_genebody_model(region_scores, regions, genes, normalize = TRUE, method.normalization = "l2")
genes       <- genes[match(rownames(gene_scores), genes$gene_id), ]
genes       <- cbind(genes, map.geneIDs(OrgDb, genes$gene_id, columns_extract = c("SYMBOL", "ENSEMBL")))

# PREPARE DATA FOR GSEA
# ---------------------
# Load the gene-set data.
cat(sprintf("Loading gene-set data from %s.\n",genesetfile))
load(genesetfile)
rownames(gene_sets) <- gene_info$Ensembl
cat(sprintf("Loaded data for %d gene sets.\n",nrow(gene_set_info)))

# Prepare the gene-set data and gene-wise statistics for the gene-set
# enrichment analysis. First, align the gene-set data with the
# gene-wise statistics.
rownames(gene_scores) <- genes$ENSEMBL

out            <- align_gene_data(gene_sets, gene_scores)
gene_sets      <- out$gene_sets
gene_scores    <- out$gene_scores
ids            <- rownames(gene_sets)
gene_info      <- gene_info[match(ids,gene_info$Ensembl),]
genes          <- genes[match(ids,genes$ENSEMBL),]
rm(out,ids)

# Next, remove gene sets with fewer than 4 genes, and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(gene_sets) >= 4 & colSums(gene_sets) <= 400)
gene_set_info <- gene_set_info[i,]
gene_sets     <- gene_sets[,i]
rm(i)

# List top 5 genes per topic
gene_scores_symbol <- gene_scores
rownames(gene_scores_symbol) <- genes[match(rownames(gene_scores), genes$ENSEMBL), "SYMBOL"]
top_genes_topics <- data.frame(matrix(nrow = 5, ncol = ncol(gene_scores_symbol)))
colnames(top_genes_topics) <- colnames(gene_scores_symbol)
for (k in colnames(gene_scores_symbol)) {
  gene_scores_topic <- gene_scores_symbol[, k]
  gene_scores_topic <- gene_scores_topic[order(abs(gene_scores_topic), decreasing = T)]
  top_genes_topics[, k] <- names(gene_scores_topic[1:5])
}

print(top_genes_topics)

# PERFORM GSEA
# ------------
# For each topic, perform a gene-set enrichment analysis using fgsea.
# Computation took 1778.35 seconds for the 13 topics.
cat("Performing gene-set enrichment analysis.\n")
timing <- system.time(
  gsea_res <- perform_gsea_all_topics(gene_sets,gene_scores,nproc = 8))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat(sprintf("Saving results to %s. \n", outfile))
save(list = c("gene_info","gene_set_info","gene_sets","genes",
              "diff_count_res","gene_scores","gene_scores_symbol","gsea_res"),
     file = outfile)
resaveRdaFiles(outfile)

# sessionInfo
sessionInfo()

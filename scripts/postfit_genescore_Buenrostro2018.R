#! /usr/bin/env Rscript
# Compute gene scores for Buenrostro et al 2018 data.
setwd("~/projects/scATACseq-topics/")
library(optparse)
library(tools)
library(Matrix)
library(fastTopics)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
source("code/gene_annotation.R")
source("code/gene_scores.R")

# Settings
# ---------
DAfile                <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=10/DA_regions_topics_noshrinkage_10000iters.rds"
genome                <- "hg19"
genescoremethod       <- "TSS"
out.dir               <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/geneanalysis-Buenrostro2018-k=10-TSS-none-l2"

cat(sprintf("DAfile               = %s \n", DAfile))
cat(sprintf("genome               = %s \n", genome))
cat(sprintf("genescoremethod      = %s \n", genescoremethod))
cat(sprintf("out.dir              = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = T)

# LOAD DATA
# ---------

# Load gene annotation
cat("Load gene annotations.\n")
if(tolower(genome) %in% c("hg19", "hg38", "mm9", "mm10")){
  cat(sprintf("load TxDb and OrgDb for %s. \n", genome))
  TxDb  <- getTxDb(genome)
  OrgDb <- getOrgDb(genome)
  genes <- get_gene_annotations(TxDb, OrgDb, columns_extract = c("ENSEMBL", "SYMBOL"))
}else{
  stop("Genome not recongized or not available. Please provide your own gene annotation data.")
}

# Prepare a data frame of gene annotation for computing gene scores,
# the first five columns need to be: chr, start, end, strand, gene_id
genes <- as.data.frame(genes)
colnames(genes)[1] <- "chr"
genes <- genes[,c("chr", "start", "end", "strand", "gene_id", "ENSEMBL", "SYMBOL")]
# Filter out genes without matching Ensembl gene IDs.
genes <- genes[!grepl("^NA_", genes$ENSEMBL), ]

# COMPUTE REGION SCORES
# ----------------------

# Load differential accessbility result.
if(!file.exists(DAfile)){
  stop("DA result not available")
}

cat("Load precomputed differential accessbility statistics.\n")
DA_res <- readRDS(DAfile)

# Filter out regions with NAs
DA_res <- DA_res[c("postmean", "z", "f0")]
rows_withNAs <- which(apply(DA_res$z, 1, anyNA))
cat("Filter out", length(rows_withNAs), "regions with NAs... \n")
DA_res$postmean <- DA_res$postmean[-rows_withNAs,]
DA_res$z <- DA_res$z[-rows_withNAs,]
DA_res$f0 <- DA_res$f0[-rows_withNAs]

# Extract genomic coordinates for ATAC-seq regions
regions <- data.frame(x = rownames(DA_res$z)) %>% tidyr::separate(x, c("chr", "start", "end"), "_") %>% dplyr::mutate_at(c("start", "end"), as.numeric)

# Compute gene-level logFC
# ---------------------------
# Compute the gene-level logFC using weighted average of region-level logFC
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene logFC using the TSS model. \n")
  gene_logFC <- compute_gene_scores_tss_model(DA_res$postmean, regions, genes, transform="none", normalization = "sum")
}else{
  cat("Compute gene logFC using the gene-body model. \n")
  gene_logFC <- compute_gene_scores_genebody_model(DA_res$postmean, regions, genes, transform="none", normalization = "sum")
}

# Compute gene scores
# ---------------------
# Compute the gene-level scores
# Combine region-level z-scores using Stouffer's z-score method
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene scores using the TSS model. \n")
  gene_scores <- compute_gene_scores_tss_model(DA_res$z, regions, genes, transform="none", normalization="l2")
}else{
  cat("Compute gene scores using the gene-body model. \n")
  gene_scores <- compute_gene_scores_genebody_model(DA_res$z, regions, genes, transform="none", normalization="l2")
}

# Compute gene-level accessbility
# -------------------------------
# Extract genomic coordinates for ATAC-seq regions
region_mean_acc <- as.matrix(DA_res$f0)

# Compute the gene-level accessbility, using weighted sum of region-level mean accessbility across topics.
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene-level mean accessbility using the TSS model. \n")
  gene_mean_acc <- compute_gene_scores_tss_model(region_mean_acc, regions, genes, transform="none", normalization = "none")[,1]
}else{
  cat("Compute gene-level mean accessbility using the gene-body model. \n")
  gene_mean_acc <- compute_gene_scores_genebody_model(region_mean_acc, regions, genes, transform="none", normalization = "none")[,1]
}

genes <- genes[match(rownames(gene_scores), genes$gene_id), ]

# Use ENSEMBL gene IDs
# rownames(gene_scores) <- genes[match(rownames(gene_scores), genes$gene_id), "ENSEMBL"]
# rownames(gene_logFC) <- genes[match(rownames(gene_logFC), genes$gene_id), "ENSEMBL"]
# names(gene_mean_acc) <- genes[match(names(gene_mean_acc), genes$gene_id), "ENSEMBL"]

genescore_res <- list(mean_acc = gene_mean_acc,
                      Z = gene_scores,
                      logFC = gene_logFC,
                      genes = genes)

saveRDS(genescore_res, file.path(out.dir, "genescore_result.rds"))

sessionInfo()

#! /usr/bin/env Rscript
# Compute gene scores based on the weighted average of region scores.
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

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--DAfile",type="character",default="DA_res.rds")
parser <- add_option(parser,"--genome",type="character",default="hg19")
parser <- add_option(parser,c("--genescoremethod","-s"),type = "character",default = "TSS")
parser <- add_option(parser,c("--transform","-t"),type = "character",default = "abs")
parser <- add_option(parser,c("--normalization","-n"),type = "character",default = "sum")
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
out    <- parse_args(parser)
DAfile               <- out$DAfile
genome               <- out$genome
genescoremethod      <- out$genescoremethod
transform.method     <- out$transform
normalization.method <- out$normalization
out.dir              <- out$out
rm(parser,out)

# ## Example settings
# DAfile            <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=11-quantile-v2/DA_regions_topics_ns1000.rds'
# genome            <- 'hg19'
# genescoremethod   <- 'TSS'
# transform.method  <- 'none'
# normalization.method  <- 'l2'
# out.dir           <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit_v2/geneanalysis-Buenrostro2018-k=11-TSS-none-l2'

cat(sprintf("DAfile               = %s \n", DAfile))
cat(sprintf("genome               = %s \n", genome))
cat(sprintf("genescoremethod      = %s \n", genescoremethod))
cat(sprintf("transform.method     = %s \n", transform.method))
cat(sprintf("normalization.method = %s \n", normalization.method))
cat(sprintf("out.dir              = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = T)

# LOAD DATA
# ---------
# Load gene annotation
# genes is a data frame containing the gene information, including: chr, start, end, strand, and gene_id.
cat("Load gene annotations.\n")
if(tolower(genome) %in% c("hg19", "hg38", "mm9", "mm10")){
  cat(sprintf("load TxDb and OrgDb for %s. \n", genome))
  TxDb  <- getTxDb(genome)
  OrgDb <- getOrgDb(genome)
  genes <- get_gene_annotations(TxDb, OrgDb, columns_extract = c("ENSEMBL", "SYMBOL"))
}else{
  stop("Genome not recongized or included. Please provide your own gene annotation data.")
}

# COMPUTE REGION SCORES USING DIFFERENTIAL ANALYSIS RESULT
# ------------------------------------------------------------

# Load differential accessbility analysis result using the topic model.
if(!file.exists(DAfile)){
  stop("DA result not available")
}

cat("Load precomputed differential accessbility statistics.\n")
DA_res <- readRDS(DAfile)

# Filter out regions with NAs
DA_res <- DA_res[c("postmean", "z", "lfsr")]
rows_withNAs <- which(apply(DA_res$z, 1, anyNA))
cat("Filter out", length(rows_withNAs), "regions with NAs... \n")
DA_res$z <- DA_res$z[-rows_withNAs,]
DA_res$postmean <- DA_res$postmean[-rows_withNAs,]
DA_res$lfsr <- DA_res$lfsr[-rows_withNAs,]

# COMPUTE GENE SCORES
# -------------------
# Prepare genes for computing gene scores, which requires the first 5 columns to be: chr, start, end, strand, gene_id
genes <- data.frame(genes)
colnames(genes)[1] <- "chr"
genes <- genes[,c("chr", "start", "end", "strand", "gene_id", "ENSEMBL", "SYMBOL")]
# Filter out genes without matching Ensembl gene ID.
genes <- genes[!grepl("^NA_", genes$ENSEMBL), ]

# Z-scores of the ATAC-seq regions from DA analysis
region_Z <- DA_res$z

# Extract genomic coordinates for ATAC-seq regions
regions <- data.frame(x = rownames(region_Z)) %>% separate(x, c("chr", "start", "end"), "_") %>% mutate_at(c("start", "end"), as.numeric)

# Compute the gene scores
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene scores using the TSS model. \n")
  gene_scores <- compute_gene_scores_tss_model(region_Z, regions, genes, transform.method=transform.method, normalization.method=normalization.method)
}else{
  cat("Compute gene scores using the gene-body model. \n")
  gene_scores <- compute_gene_scores_genebody_model(region_Z, regions, genes, transform.method=transform.method, normalization.method=normalization.method)
}

genes <- genes[match(rownames(gene_scores), genes$gene_id), ]

# COMPUTE GENE LEVEL logFC
# ---------------------------
region_logFC <- DA_res$postmean

# Compute the gene logFC
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene logFC using the TSS model. \n")
  gene_logFC <- compute_gene_scores_tss_model(region_logFC, regions, genes, transform.method=transform.method, normalization.method = "sum")
}else{
  cat("Compute gene logFC using the gene-body model. \n")
  gene_logFC <- compute_gene_scores_genebody_model(region_logFC, regions, genes, transform.method=transform.method, normalization.method = "sum")
}

if(!all.equal(rownames(gene_logFC), genes$gene_id))
  stop("ERROR: Gene names do not match!")


# rownames(gene_scores) <- genes[match(rownames(gene_scores), genes$gene_id), "ENSEMBL"]
# rownames(gene_logFC) <- genes[match(rownames(gene_logFC), genes$gene_id), "ENSEMBL"]

genescore_res <- list(Z = gene_scores,
                      logFC = gene_logFC,
                      genes = genes)

saveRDS(genescore_res, file.path(out.dir, "genescore_result.rds"))

# sessionInfo()

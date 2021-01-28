#! /usr/bin/env Rscript
# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# compute gene scores based on the weighted average of region scores,
# and perform gene-set enrichment analysis based on a multinomial topic model.

library(optparse)
library(tools)
library(Matrix)
library(fastTopics)
library(fgsea)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
source("~/projects/scATACseq-topics/code/gene_annotation.R")
source("~/projects/scATACseq-topics/code/gene_scores.R")
source("~/projects/scATACseq-topics/code/gsea.R")

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--modelfit",type = "character",default="fit.rds")
parser <- add_option(parser,"--geneset",type="character",default="gene_sets.RData")
parser <- add_option(parser,"--genome",type="character",default="hg19")
parser <- add_option(parser,c("--genescoremethod","-s"),type = "character",default = "genebody")
parser <- add_option(parser,c("--normalization","-n"),type = "character",default = "sum")
parser <- add_option(parser,"--runGSEA",action = "store_true",default = FALSE)
parser <- add_option(parser,"--nc",type = "integer",default = 1)
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
out    <- parse_args(parser)
countsfile      <- out$counts
modelfitfile    <- out$modelfit
genesetfile     <- out$geneset
genome          <- out$genome
genescoremethod <- out$genescoremethod
normalization   <- out$normalization
runGSEA         <- out$runGSEA
nc              <- out$nc
out.dir         <- out$out
rm(parser,out)

cat(sprintf("countsfile      = %s \n", countsfile))
cat(sprintf("modelfitfile    = %s \n", modelfitfile))
cat(sprintf("genesetfile     = %s \n", genesetfile))
cat(sprintf("genome          = %s \n", genome))
cat(sprintf("genescoremethod = %s \n", genescoremethod))
cat(sprintf("normalization   = %s \n", normalization))
cat(sprintf("runGSEA         = %s \n", runGSEA))
cat(sprintf("nc              = %s \n", nc))
cat(sprintf("out.dir         = %s \n", out.dir))

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

# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# LOAD MODEL FIT
# --------------
cat(sprintf("Loading Poisson NMF model fit from %s\n",modelfitfile))
fit <- readRDS(modelfitfile)$fit

# COMPUTE REGION SCORES USING DIFFERENTIAL ANALYSIS
# -------------------------------------------------
# Perform differential accessbility analysis using the multinomial topic model.
outfile <- file.path(out.dir, "diffcount_regions_topics.rds")
if(file.exists(outfile)){
  cat("Load precomputed differential accessbility statistics.\n")
  diff_count_res <- readRDS(outfile)
}else{
  cat("Computing differential accessbility statistics from topic model.\n")
  timing <- system.time(diff_count_res <- diff_count_analysis(fit,counts))
  cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
  cat("Saving results.\n")
  saveRDS(diff_count_res, outfile)
}

# COMPUTE GENE SCORES
# -------------------
# Prepare genes for computing gene scores, which requires the first 5 columns to be: chr, start, end, strand, gene_id
genes <- data.frame(genes)
colnames(genes)[1] <- "chr"
genes <- genes[,c("chr", "start", "end", "strand", "gene_id", "ENSEMBL", "SYMBOL")]
# Filter out genes without matching Ensembl gene ID.
genes <- genes[!grepl("^NA_", genes$ENSEMBL), ]

# Extract genomic coordinates for ATAC-seq regions
region_Z <- diff_count_res$Z
regions <- data.frame(x = rownames(region_Z)) %>% separate(x, c("chr", "start", "end"), "_") %>% mutate_at(c("start", "end"), as.numeric)

# Compute the gene scores (z-scores)
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene z-scores using the TSS model. \n")
  gene_scores <- compute_gene_scores_tss_model(region_Z, regions, genes, normalize = TRUE, method.normalization = normalization)
}else{
  cat("Compute gene z-scores using the gene-body model. \n")
  gene_scores <- compute_gene_scores_genebody_model(region_Z, regions, genes, normalize = TRUE, method.normalization = normalization)
}

genes <- genes[match(rownames(gene_scores), genes$gene_id), ]

# COMPUTE GENE weighted logFC
# ---------------------------
# Extract genomic coordinates for ATAC-seq regions
region_beta <- diff_count_res$beta

# Compute the gene logFC
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene logFC using the TSS model. \n")
  gene_logFC <- compute_gene_scores_tss_model(region_beta, regions, genes, normalize = TRUE, method.normalization = "sum")
}else{
  cat("Compute gene logFC using the gene-body model. \n")
  gene_logFC <- compute_gene_scores_genebody_model(region_beta, regions, genes, normalize = TRUE, method.normalization = "sum")
}

if(!all.equal(rownames(gene_logFC), genes$gene_id))
  stop("ERROR: Gene names do not match!")

# COMPUTE GENE weighted average accessbility
# -----------------------------------------------------
# Extract genomic coordinates for ATAC-seq regions
region_mean <- as.matrix(diff_count_res$colmeans)

# Compute the gene mean accessbility, by weighted sum of the mean accessbility across topics.
if(toupper(genescoremethod) == "TSS"){
  cat("Compute gene mean accessbility using the TSS model. \n")
  gene_mean_acc <- compute_gene_scores_tss_model(region_mean, regions, genes, normalize = TRUE, method.normalization = "sum")[,1]
}else{
  cat("Compute gene mean accessbility using the gene-body model. \n")
  gene_mean_acc <- compute_gene_scores_genebody_model(region_mean, regions, genes, normalize = TRUE, method.normalization = "sum")[,1]
}

if(!all.equal(names(gene_mean_acc), genes$gene_id))
  stop("ERROR: Gene names do not match!")

# rownames(gene_scores) <- genes[match(rownames(gene_scores), genes$gene_id), "ENSEMBL"]
# rownames(gene_logFC) <- genes[match(rownames(gene_logFC), genes$gene_id), "ENSEMBL"]
# names(gene_mean_acc) <- genes[match(names(gene_mean_acc), genes$gene_id), "ENSEMBL"]

genescore_res <- list(colmeans = gene_mean_acc,
                      Z = gene_scores,
                      beta = gene_logFC,
                      genes = genes)

saveRDS(genescore_res, file.path(out.dir, "genescore_result_topics.rds"))

# GENE SET ENRICHMENT ANALYSIS
# ----------------------------
if(runGSEA) {

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
  rownames(gene_scores) <- genes[match(rownames(gene_scores), genes$gene_id), "ENSEMBL"]

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

  # PERFORM GSEA
  # ------------
  # For each topic, perform a gene-set enrichment analysis using fgsea.
  cat("Performing gene-set enrichment analysis.\n")
  timing <- system.time(
    gsea_res <- perform_gsea_all_topics(gene_sets,gene_scores,nproc = nc))
  cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

  outfile <- file.path(out.dir, "genescores_gsea.Rdata")
  cat(sprintf("Saving gene scores and GSEA results to %s. \n", outfile))
  save(list = c("gene_info","gene_set_info","gene_sets",
                "genes","gene_scores","gsea_res"),
       file = outfile)
  resaveRdaFiles(outfile)

}

# sessionInfo
sessionInfo()

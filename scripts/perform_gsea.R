#! /usr/bin/env Rscript
# Perform gene-set enrichment analysis based on gene scores from topic model.

library(optparse)
library(tools)
library(Matrix)
library(fastTopics)
library(pathways)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--genescore",type="character",default=NA)
parser <- add_option(parser,"--genome",type="character",default=NA)
parser <- add_option(parser,"--geneset",type="character",default=NA)
parser <- add_option(parser,"--nc",type="integer",default = 1)
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
out    <- parse_args(parser)
genescorefile  <- out$genescore
genome         <- out$genome
genesetfile    <- out$geneset
nc             <- out$nc
out.dir        <- out$out
rm(parser,out)

# genescorefile <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit/geneanalysis-Buenrostro2018-k=11-TSS-none-l2/genescore_result_topics.rds"
# genome <- "hg19"
# out.dir <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit/geneanalysis-Buenrostro2018-k=11-TSS-none-l2"

cat(sprintf("genescorefile = %s \n", genescorefile))
cat(sprintf("genome        = %s \n", genome))
cat(sprintf("nc            = %s \n", nc))
cat(sprintf("out.dir       = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = T)

# LOAD DATA
# ---------

# LOAD GENE SCORES
# -------------------
genescore_res <- readRDS(genescorefile)
genes <- genescore_res$genes
gene_scores <- genescore_res$Z
rownames(gene_scores) <- genes[match(rownames(gene_scores), genes$gene_id), "ENSEMBL"]

# GENE SET ENRICHMENT ANALYSIS
# ----------------------------

# PREPARE DATA FOR GSEA
# ---------------------
# Load the gene-set data.
if (genome %in% c("hg19", "hg38", "human")) {
  cat("Loading human gene set data.\n")
  data(gene_sets_human)
  gene_sets <- gene_sets_human$gene_sets
  gene_set_info <- gene_sets_human$gene_set_info
} else if (genome %in% c("mm9", "mm10", "mouse")) {
  cat("Loading mouse gene set data.\n")
  data(gene_sets_mouse)
  gene_sets <- gene_sets_mouse$gene_sets
  gene_set_info <- gene_sets_mouse$gene_set_info
} else{
  cat(sprintf("Loading gene-set data from %s.\n",genesetfile))
  load(genesetfile)
}

# Remove gene sets with fewer than 4 genes, and with more than
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
  gsea_res <- pathways::perform_gsea(gene_sets, gene_scores, nproc = nc))

cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

outfile <- file.path(out.dir, "gsea_result.rds")
cat(sprintf("Saving GSEA results to %s \n", outfile))
saveRDS(gsea_res, outfile)

# PLOT GSEA RESULT
# ----------------

# Create interactive plots for exploring the GSEA results.
cat("Making interactive GSEA plots... \n")

for ( k in colnames(gsea_res$pval)) {
  outfile <- paste0(out.dir, "/gsea_plots/gsea_plotly_", k, ".html")
  gsea_plotly(gsea_res,gene_set_info, k, file = outfile, title = k)
}

# sessionInfo
sessionInfo()

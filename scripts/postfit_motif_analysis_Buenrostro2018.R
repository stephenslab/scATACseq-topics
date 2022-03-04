#! /usr/bin/env Rscript
# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis using HOMER.

setwd("~/projects/scATACseq-topics/")
library(optparse)
library(Matrix)
library(fastTopics)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
source("code/motif_analysis.R")

## Settings
DAfile          <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=10/DA_regions_topics_noshrinkage_10000iters.rds'
genome          <- 'hg19'
homerpath       <- '/project2/xinhe/software/homer/bin/findMotifsGenome.pl'
nc              <- 8
out.dir         <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/motifanalysis-Buenrostro2018-k=10'

cat(sprintf("DAfile       = %s \n", DAfile))
cat(sprintf("genome       = %s \n", genome))
cat(sprintf("homerpath    = %s \n", homerpath))
cat(sprintf("nc           = %s \n", nc))
cat(sprintf("out.dir      = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# LOAD DATA
# ---------

# DIFFERENTIAL ACCESSBILITY ANALYSIS
# ------------------------------------------
# Load differential accessbility analysis result using the topic model.
if(!file.exists(DAfile)){
  stop("DA result not available")
}

cat("Load precomputed differential accessbility statistics.\n")
DA_res <- readRDS(DAfile)

# Regions with NAs
rows_withNAs <- which(apply(DA_res$z, 1, anyNA))
cat(length(rows_withNAs), "regions with NAs in z-scores... \n")
head(DA_res$z[rows_withNAs,])

# Filter out regions with NAs in z-scores
# DA_res$z <- DA_res$z[-rows_withNAs,]
# DA_res$postmean <- DA_res$postmean[-rows_withNAs,]
# DA_res$lpval <- DA_res$lpval[-rows_withNAs,]
# DA_res$lfsr <- DA_res$lfsr[-rows_withNAs,]

# SELECT REGIONS and PERFORM MOTIF ENRICHMENT ANALYSIS USING HOMER
# -----------------------------------------------------------------

# Select regions with DA p-value < 0.1
# -------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_pval_0.1_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "pval", thresh.pval = 0.1, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

# For each topic, perform TF motif enrichment analysis using HOMER hypergeometric test.
cat("Performing motif enrichment analysis using HOMER.\n")
homer_res <- vector("list", ncol(DA_res$z))
names(homer_res) <- colnames(DA_res$z)
for(k in 1:ncol(DA_res$z)){
  homer_res[[k]] <- run_homer(selected_regions$filenames[k],
                              genome = genome,
                              homer.path = homerpath,
                              use.hypergeometric = TRUE,
                              out.dir=paste0(homer.dir, "/homer_result_topic_", k),
                              n.cores=nc)
}
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))


# Select regions with DA p-value < 0.05
# --------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_pval_0.05_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "pval", thresh.pval = 0.05, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

# For each topic, perform TF motif enrichment analysis using HOMER hypergeometric test.
cat("Performing motif enrichment analysis using HOMER.\n")
homer_res <- vector("list", ncol(DA_res$z))
names(homer_res) <- colnames(DA_res$z)
for(k in 1:ncol(DA_res$z)){
  homer_res[[k]] <- run_homer(selected_regions$filenames[k],
                              genome = genome,
                              homer.path = homerpath,
                              use.hypergeometric = TRUE,
                              out.dir=paste0(homer.dir, "/homer_result_topic_", k),
                              n.cores=nc)
}
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))


sessionInfo()

#! /usr/bin/env Rscript
# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis based on a multinomial topic model.

library(optparse)
library(Matrix)
library(fastTopics)
library(dplyr)
library(tidyr)
source("~/projects/scATACseq-topics/code/motif_analysis.R")

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--modelfit",type = "character",default="fit.rds")
parser <- add_option(parser,"--genome",type="character",default="hg19")
parser <- add_option(parser,"--homerpath",type="character",default="findMotifsGenome.pl")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
parser <- add_option(parser,"--testrun", action="store_true", default = FALSE)
out    <- parse_args(parser)
countsfile      <- out$counts
modelfitfile    <- out$modelfit
genome          <- out$genome
homerpath       <- out$homerpath
nc              <- out$nc
out.dir         <- out$out
testrun         <- out$testrun
rm(parser,out)

if(testrun){
  countsfile   <- "/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data_Chen2019pipeline/Buenrostro_2018_binarized_counts.RData"
  modelfitfile <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/fit-Buenrostro2018-binarized-scd-ex-k=11.rds"
  genome       <- "hg19"
  homerpath    <- "/project2/xinhe/software/homer/bin/findMotifsGenome.pl"
  nc           <- 8
  out.dir      <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized"
}

cat(sprintf("countsfile   = %s \n", countsfile))
cat(sprintf("modelfitfile = %s \n", modelfitfile))
cat(sprintf("genome       = %s \n", genome))
cat(sprintf("homerpath    = %s \n", homerpath))
cat(sprintf("nc           = %s \n", nc))
cat(sprintf("out.dir      = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = T)

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

# PERFORM DIFFERENTIAL ACCESSBILITY ANALYSIS
# ------------------------------------------
# Perform differential accessbility analysis using the multinomial topic model.
cat("Computing differential accessbility statistics from topic model.\n")
outfile <- file.path(out.dir, "diffcount_regions_topics.rds")
if(file.exists(outfile)){
  diff_count_res <- readRDS(outfile)
}else{
  timing <- system.time(diff_count_res <- diff_count_analysis(fit,counts))
  cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
  cat("Saving results.\n")
  saveRDS(diff_count_res, outfile)
}

# SELECT REGIONS FOR MOTIF ENRICHMENT ANALYSIS
# ------------------------------------------
cat("Select regions. \n")
homer.dir <- paste0(out.dir, "/HOMER")
selected_regions <- select_regions(diff_count_res, method="FDR&logFC", out.dir = homer.dir,
                           thresh.FDR = 0.05, thresh.logFC = 1, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions_FDR0.05_logFC1.rds"))

# PERFORM MOTIF ENRICHMENT ANALYSIS USING HOMER
# ---------------------------------------------
# For each topic, perform a gene-set enrichment analysis using fgsea.
cat("Performing motif enrichment analysis.\n")
homer_res <- vector("list", ncol(diff_count_res$Z))
names(homer_res) <- colnames(diff_count_res$Z)
for(k in 1:ncol(diff_count_res$Z)){
  homer_res[[k]] <- run_homer(selected_regions$filenames[k],
                              genome,
                              homerpath,
                              out.dir=paste0(homer.dir, "/homer_result_topic_", k),
                              n.cores=nc)
}
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults_topics.rds"))

# sessionInfo
sessionInfo()

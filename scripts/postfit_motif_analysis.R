#! /usr/bin/env Rscript
# Perform differential accessbility analysis for ATAC-seq regions (peaks),
# and perform TF motif enrichment analysis based on a multinomial topic model.

setwd("~/projects/scATACseq-topics/")
library(optparse)
library(Matrix)
library(fastTopics)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
source("code/motif_analysis.R")

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--DAfile",type="character",default="DA_res.rds")
parser <- add_option(parser,"--genome",type="character",default="hg19")
parser <- add_option(parser,"--selectmethod",type="character",default="lfsr")
parser <- add_option(parser,"--homerpath",type="character",default="findMotifsGenome.pl")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
out    <- parse_args(parser)
DAfile          <- out$DAfile
genome          <- out$genome
selectmethod    <- out$selectmethod
homerpath       <- out$homerpath
nc              <- out$nc
out.dir         <- out$out
rm(parser,out)

## Example settings
# DAfile       <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=11-quantile-v2/DA_regions_topics_ns1000.rds'
# genome          <- 'hg19'
# selectmethod    <- 'lfsr'
# homerpath       <- '/project2/xinhe/software/homer/bin/findMotifsGenome.pl'
# nc              <- 8
# out.dir         <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018_Chen2019pipeline/binarized/postfit_v2/motifanalysis-Buenrostro2018-k=11-quantile'

cat(sprintf("DAfile    = %s \n", DAfile))
cat(sprintf("genome       = %s \n", genome))
cat(sprintf("homerpath    = %s \n", homerpath))
cat(sprintf("selectmethod = %s \n", selectmethod))
cat(sprintf("nc           = %s \n", nc))
cat(sprintf("out.dir      = %s \n", out.dir))

if(!dir.exists(out.dir))
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

# LOAD DATA
# ---------

# PERFORM DIFFERENTIAL ACCESSBILITY ANALYSIS
# ------------------------------------------
# Load differential accessbility analysis result using the topic model.
if(!file.exists(DAfile)){
  stop("DA result not available")
}

cat("Load precomputed differential accessbility statistics.\n")
DA_res <- readRDS(DAfile)

# SELECT REGIONS FOR MOTIF ENRICHMENT ANALYSIS
# --------------------------------------------
cat("Select regions. \n")
homer.dir <- paste0(out.dir, "/HOMER/")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = 'quantile', thresh.quantile = 0.99, out.dir = homer.dir, save.bed = TRUE)

saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

# PERFORM MOTIF ENRICHMENT ANALYSIS USING HOMER
# ---------------------------------------------
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

# sessionInfo()

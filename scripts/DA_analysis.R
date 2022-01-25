#! /usr/bin/env Rscript
# Perform differential accessbility (DA) analysis for ATAC-seq regions,
# compute gene scores based on the weighted average of region scores.

library(optparse)
library(Matrix)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--modelfit",type = "character",default="fit.rds")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
parser <- add_option(parser,"--ns",type = "integer",default = 1000)
parser <- add_option(parser,"--nsplit",type = "integer",default = 100)
parser <- add_option(parser,c("--out","-o"),type="character",default="out")
out    <- parse_args(parser)

countsfile           <- out$counts
modelfitfile         <- out$modelfit
nc                   <- out$nc
ns                   <- out$ns
nsplit               <- out$nsplit
outdir               <- out$out
rm(parser,out)

cat(sprintf("countsfile           = %s \n", countsfile))
cat(sprintf("modelfitfile         = %s \n", modelfitfile))
cat(sprintf("nc                   = %s \n", nc))
cat(sprintf("ns                   = %s \n", ns))
cat(sprintf("nsplit               = %s \n", nsplit))
cat(sprintf("outdir              = %s \n", outdir))

if(!dir.exists(outdir))
  dir.create(outdir, showWarnings = FALSE, recursive = T)

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

# COMPUTE REGION SCORES USING DIFFERENTIAL ANALYSIS
# -------------------------------------------------
# Perform differential accessbility analysis using the multinomial topic model.
outfile <- file.path(outdir, paste0("DA_regions_topics_", ns,"iters.rds"))

if(file.exists(outfile)){
  cat("Load precomputed differential accessbility statistics.\n")
  DA_res <- readRDS(outfile)
}else{
  cat("Computing differential accessbility statistics from topic model.\n")
  cat("Run", ns, "iterations of MCMC...\n")
  timing <- system.time(
    DA_res <- de_analysis(fit,counts,pseudocount = 0.1,
                        control = list(ns = ns,nc = nc,nsplit = nsplit)))
  cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
  cat("Saving results.\n")
  saveRDS(DA_res, outfile)
}

# sessionInfo()

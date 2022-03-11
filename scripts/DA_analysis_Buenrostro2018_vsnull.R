#! /usr/bin/env Rscript
# Perform differential accessbility (DA) analysis for Buenrostro 2018 scATAC-seq data based on topic modeling fit

library(optparse)
library(Matrix)
library(fastTopics)

countsfile           <- "/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data/Buenrostro_2018_binarized.RData"
modelfitfile         <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/fit-Buenrostro2018-binarized-scd-ex-k=10.rds"
nc                   <- 20
ns                   <- 10000
nsplit               <- 100
outdir               <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=10"

cat(sprintf("countsfile   = %s \n", countsfile))
cat(sprintf("modelfitfile = %s \n", modelfitfile))
cat(sprintf("nc           = %s \n", nc))
cat(sprintf("ns           = %s \n", ns))
cat(sprintf("nsplit       = %s \n", nsplit))
cat(sprintf("outdir       = %s \n", outdir))

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

# DIFFERENTIAL ACCESSIBILITY ANALYSIS
# -------------------------------------------------

# Perform differential accessbility for logFC relative to the null (average) with ash shrinkage
outfile <- file.path(outdir, paste0("DA_regions_topics_vsnull_ash_", ns,"iters.rds"))

cat("Computing differential accessbility statistics from topic model.\n")
cat("Run", ns, "iterations of MCMC...\n")
timing <- system.time(
  DA_res <- de_analysis(fit,counts,pseudocount = 0.1, lfc.stat = "vsnull",
                        control = list(ns = ns,nc = nc,nsplit = nsplit)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("Saving results.\n")
saveRDS(DA_res, outfile)

# Perform differential accessbility for logFC relative to the null (average) without shrinkage
outfile <- file.path(outdir, paste0("DA_regions_topics_vsnull_noshrinkage_", ns,"iters.rds"))

cat("Computing differential accessbility statistics from topic model.\n")
cat("Run", ns, "iterations of MCMC...\n")
timing <- system.time(
  DA_res <- de_analysis(fit,counts,pseudocount = 0.1,shrink.method = "none",
                        lfc.stat = "vsnull",
                        control = list(ns = ns,nc = nc,nsplit = nsplit)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("Saving results.\n")
saveRDS(DA_res, outfile)

sessionInfo()

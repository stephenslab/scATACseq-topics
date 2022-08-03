#! /usr/bin/env Rscript
#
# Fit a Poisson non-negative factorization to the Cusanovich et al (2018) single-cell ATAC-seq data set.

# Load packages.
library(optparse)
library(Matrix)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--tissue",type="character",default="")
parser <- add_option(parser,c("--k","-k"),type = "integer",default = 10)
parser <- add_option(parser,c("--numprefit"),type="integer",default=200)
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=300)
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
# print(out)
countsfile  <- out$counts
tissue      <- out$tissue
k           <- out$k
numprefit   <- out$numprefit
numiter     <- out$numiter
nc          <- out$nc
rm(parser,out)

method <- "scd"
extrapolate <- TRUE
outdir <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Cusanovich_2018/tissues"
prefitfile <- file.path(outdir, sprintf("prefit-Cusanovich2018-%s-scd-ex-k=%d.rds", tissue, k))
outfile <- file.path(outdir, sprintf("fit-Cusanovich2018-%s-scd-ex-k=%d.rds", tissue, k))

cat("prefitfile:", prefitfile, "\n")
cat("outfile:", outfile, "\n")

if (!dir.exists(outdir))
  dir.create(outdir)

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

if(any(colSums(counts) < 20)){
  # Filtered out peaks with accessbility in fewer than 20 cells
  cat("Filter out peaks with accessbility in fewer than 20 cells...\n")
  i      <- which(colSums(counts) >= 20)
  counts <- counts[,i]
  cat(sprintf("After filtering, we have %d rows and %d columns left. \n",nrow(counts),ncol(counts)))
  save(list = c("samples","counts"), file = countsfile)
}

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# PREFIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ------------------------------------------------
cat(sprintf("Running %d EM updates to identify a good initialization.\n",
            numprefit))
timing <- system.time(
  fit0 <- fit_poisson_nmf(counts,k = k,numiter = numprefit,method = "em",
                         control = list(numiter = 4,nc = nc)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE PREFIT RESULTS
# -------------------
cat("Saving prefit results.\n")
saveRDS(list(k = k,fit = fit0),file = prefitfile)

# FIT POISSON NMF MODEL
# ---------------------------------------------
cat("Fitting Poisson NMF to count data.\n")
control <- list(extrapolate = TRUE,nc = nc,
                numiter = ifelse(method == "ccd",1,4))
timing <- system.time({
  fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = numiter,
                         method = method,control = control)
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat("Saving Poisson NMF fit results.\n")
saveRDS(list(method = method,control = control,fit = fit),
        file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())

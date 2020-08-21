#! /usr/bin/env Rscript
#
# Script for "prefitting" a Poisson non-negative factorization to one
# of the single-cell RNA-seq data sets.
#
# This script is intended to be run from the command-line shell, with
# options that are processed using the optparse package. For example,
# to fit a rank-4 Poisson non-negative matrix factorization to counts
# data test.RData by running 500 EM updates, 2 threads, and with
# results saved to prefit_test.rds, run this command:
#
#   ./prefit_poisson_nmf_purified_pbmc.R --counts test.RData \
#      -k 4 -n 500 --nc 2 -o prefit_test.rds
#
# Running the script without specifying any options will prefit a
# rank-3 Poisson non-negative matrix factorization to data set
# counts.RData by running 1,000 EM updates without multithreading, and
# with results saved to out.rds.
#
# The input .RData file specified by --counts should contain a matrix,
# "counts", containing the count data that will be provided as input
# to fit_poisson_nmf.
#

# Load a few packages.
library(Matrix)
library(optparse)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,c("--k","-k"),type = "integer",default = 3)
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=1000)
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countsfile  <- out$counts
outfile     <- out$out
k           <- out$k
numiter     <- out$numiter
nc          <- out$nc
rm(parser,out)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# PREFIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ------------------------------------------------
# The aim here is to run enough EM updates so that it is difficult for
# the other algorithms to "escape" this local maximum of the
# likelihood surface.
cat(sprintf("Running %d EM updates to identify a good initialization.\n",
            numiter))
timing <- system.time(
  fit <- fit_poisson_nmf(counts,k = k,numiter = numiter,method = "em",
                         control = list(numiter = 4,nc = nc)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(list(k = k,fit = fit),file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())

#! /usr/bin/env Rscript
#
# Fit a Poisson non-negative factorization to one of the single-cell
# RNA-seq data sets.
#
# This script is intended to be run from the command-line shell, with
# options that are processed with the optparse package. For example,
# to fit a Poisson non-negative matrix factorization to counts data
# test.RData by running 500 SCD updates with extrapolation, 2 threads,
# with estimates initialized from test_prefit.rds, and with results
# saved to test_fit.rds, run this command:
#
#   ./fit_poisson_nmf.R --counts test.RData --prefit test_prefit.rds \
#     --method scd -n 500 --extrapolate --nc 2 -o test_fit.rds
#
# Running the script without specifying any options will fit a Poisson
# non-negative matrix factorization by running 1,000 EM updates
# without extrapolation, without multithreading, initialized to
# prefit.rds, and with results saved to out.rds.
#
# The input .RData file specified by --counts should contain a matrix,
# "counts", containing the count data that will be provided as input
# to fit_poisson_nmf. The input .rds file specified by --prefit should
# contain an object "fit", an output from a previous call to
# fit_poisson_nmf or init_poisson_nmf.
#

# Load a few packages.
library(Matrix)
library(optparse)
library(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,"--prefit",type = "character",default="prefit.rds")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,"--method",type = "character",default = "em")
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=1000)
parser <- add_option(parser,"--extrapolate",action = "store_true")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countsfile  <- out$counts
prefitfile  <- out$prefit
outfile     <- out$out
method      <- out$method
numiter     <- out$numiter
extrapolate <- !is.null(out$extrapolate)
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

# LOAD PRE-FITTED MODEL
# ---------------------
cat(sprintf("Loading pre-fitted model from %s\n",prefitfile))
fit0 <- readRDS(prefitfile)$fit

# FIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ---------------------------------------------
# Now we are ready to perform the main model-fitting step.
cat("Fitting Poisson NMF to count data.\n")
control <- list(extrapolate = extrapolate,nc = nc,
                numiter = ifelse(method == "ccd",1,4))
timing <- system.time({
  fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = numiter,
                         method = method,control = control)
})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(list(method = method,control = control,fit = fit),
        file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())

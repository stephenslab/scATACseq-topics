#! /usr/bin/env Rscript
#
# Fit a Poisson non-negative factorization to the Buenrostro et al
# (2018) single-cell ATAC-seq data set.
# These were the steps taken to
# load R and allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 10 --mem=20G --time=24:00:00
#   module load R/4.0.4
#
#

# Load packages.
library(Matrix)
library(fastTopics)

# Process the command-line arguments.
countsfile <- "/project2/mstephens/kevinluo/scATACseq-topics/data/Buenrostro_2018/processed_data/Buenrostro_2018_binarized.RData"
method <- "scd"
k <- 8
nc <- 10
num.prefit.iter <- 200
numiter <- 500
outdir <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/filtered_peaks"
prefitfile <- file.path(outdir, "prefit-Buenrostro2018-binarized-k=8")
outfile <- file.path(outdir, "fit-Buenrostro2018-binarized-k=8")

if (!dir.exists(outdir))
  dir.create(outdir)

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# Filter out chromatin accessibility regions with low mean accessibility.
i      <- which(colSums(counts) >= 20)
peaks  <- peaks[i,]
counts <- counts[,i]
cat(sprintf("After filtering, we have %d rows and %d columns left. \n",nrow(counts),ncol(counts)))
save(list = c("samples", "peaks", "counts"),
     file = file.path(outdir, "Buenrostro_2018_binarized_filtered.RData"))

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# PREFIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ------------------------------------------------
cat(sprintf("Running %d EM updates to identify a good initialization.\n",
            num.prefit.iter))
timing <- system.time(
  fit0 <- fit_poisson_nmf(counts,k = k,numiter = num.prefit.iter,method = "em",
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

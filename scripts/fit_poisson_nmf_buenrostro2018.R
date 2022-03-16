# Fit a Poisson non-negative factorization to the Buenrostro et al
# (2018) single-cell ATAC-seq data set. These were the steps taken to
# load R and allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 8 --mem=16G --time=24:00:00
#   module load R/3.5.1
# 
library(Matrix)
library(fastTopics)

# Load the previously prepared count data.
load("../data/Buenrostro_2018_binarized.RData")

# Filter out chromatin accessibility regions with low mean accessibility.
i      <- which(colSums(counts) >= 20)
peaks  <- peaks[i,]
counts <- counts[,i]

# Fit a Poisson NMF model.
t0 <- proc.time()
fit <- fit_poisson_nmf(counts,k = 8,numiter = 200,method = "scd",
                       init.method = "random",verbose = "detailed",
                       control = list(numiter = 4,nc = 8,extrapolate = TRUE))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
saveRDS(fit,file = "fit-buenrostro-2018-k=8.rds")

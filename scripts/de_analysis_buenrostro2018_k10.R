# Perform GoM DE analysis using the binomial topic model fitted to the
# Buenrostro et al (2018) data, with k = 10 topics. These were the
# steps taken to load R and allocate computing resources for this
# analysis:
#
#   sinteractive -p broadwl -c 8 --mem=16G --time=36:00:00
#   module load R/3.5.1
#

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data.
load(file.path("../output/Buenrostro_2018/binarized/filtered_peaks",
               "Buenrostro_2018_binarized_filtered.RData"))

# Load the k = 10 binomial topic model fit.
load(file.path("../output/Buenrostro_2018/binarized/filtered_peaks",
               "binom-fit-buenrostro2018-k=10.RData"))

# Perform the DE analysis.
t0 <- proc.time()
de <- de_analysis(fit_binom,counts,shrink.method = "none",pseudocount = 0.1,
                  control = list(ns = 1e5,nc = 8))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = "de",file = "de-buenrostro2018-k=10-noshrink.RData")
resaveRdaFiles("de-buenrostro2018-k=10-noshrink.RData")

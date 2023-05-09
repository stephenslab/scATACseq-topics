# TO DO: Explain here what this script is for, and how to use it.
#
#   sinteractive -p broadwl -c 8 --mem=32G --time=10:00:00
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

# Load the k = 10 Poisson NMF model fit.
fit_pois <- readRDS(
  file.path("../output/Buenrostro_2018/binarized/filtered_peaks",
            "fit-Buenrostro2018-binarized-filtered-scd-ex-k=10.rds"))$fit

# Perform the conversion a second time, this time with several EM
# updates to refine the fit.
t0 <- proc.time()
fit_binom_em <- poisson2binom(counts,fit_pois,numem = 4)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("fit_binom","fit_binom_em"),
     file = "binom-fit-buenrostro2018-k=10.RData")
resaveRdaFiles("binom-fit-buenrostro2018-k=10.RData")


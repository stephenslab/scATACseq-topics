# TO DO: Explain here what this script is for, and how to use it.

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
fit <- readRDS(
  file.path("../output/Buenrostro_2018/binarized/filtered_peaks",
            "fit-Buenrostro2018-binarized-filtered-scd-ex-k=10.rds"))$fit
fit <- poisson2multinom(fit)

# Perform the DE analysis.
t0 <- proc.time()
timing <- system.time(
  DA_res <- de_analysis(fit,counts,shrink.method = "none",pseudocount = 0.1,
                        control = list(ns = 1000,nc = 4)))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = "de",file = "de-buenrostro2018-k=10-noshrink.RData")
resaveRdaFiles("de-buenrostro2018-k=10-noshrink.RData")

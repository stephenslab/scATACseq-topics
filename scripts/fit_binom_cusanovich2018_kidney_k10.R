# TO DO: Explain here what this script is for, and how to use it.

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data.
load(file.path("../data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))

# Load the k = 10 Poisson NMF model fit.
fit <- readRDS(
  file.path("../output/Cusanovich_2018/tissues",
            "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit

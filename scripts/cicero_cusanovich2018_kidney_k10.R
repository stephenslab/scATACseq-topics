# TO DO: Explain here what this script does, and how to use it.

# Load a few packages.
library(tools)
library(Matrix)
library(cicero)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data, and reformat it as a
# "CellDataSet" (CDS) object.
load(file.path("../data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))
counts <- t(counts)
dat    <- summary(counts)
dat$i  <- rownames(counts)[dat$i]
dat$j  <- colnames(counts)[dat$j]
names(dat) <- c("peak","cell","count")
counts <- dat
rm(dat)
stop()
cds <- make_atac_cds(counts,binarize = TRUE)

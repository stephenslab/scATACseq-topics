#! /usr/bin/env Rscript

## Fit cisTopic

# Load packages
library(optparse)
library(Matrix)
suppressPackageStartupMessages(library(cisTopic))

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--counts",type="character",default="counts.RData")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,c("--k","-k"),type = "integer",default = 3)
parser <- add_option(parser,c("--numiter","-n"),type="integer",default=500)
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countsfile  <- out$counts
outfile     <- out$out
k           <- out$k
numiter     <- out$numiter
nc          <- out$nc
rm(parser,out)

cat(sprintf("cisTopic analysis using data: %s, with %d topics \n", countsfile, k))

# Initialize the sequence of pseudorandom numbers.
seed <- 1

# LOAD DATA
# ---------
# Load the previously prepared count data.
cat(sprintf("Loading data from %s.\n",countsfile))
load(countsfile)
cat(sprintf("Loaded %d x %d counts matrix.\n",nrow(counts),ncol(counts)))

# CREATE cisTopicObject
# ---------------------------------------------
# The precomputed matrix should have cells as columns, regions as rows and fragments/reads counts as values.
counts <- t(counts)
# The rownames of the matrix must contain the region coordinates in position format (e.g. chr1:123456-134567)
rownames(counts) <- paste0(peaks$chr, ":", peaks$start, "-",peaks$end)

# Initialize the cisTopic object
project.name <- tools::file_path_sans_ext(basename(outfile))
cisTopicObject <- createcisTopicObject(counts, project.name = project.name)

# add cell metadata
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = samples)

# FIT cisTopic models
# -------------------
# Fit cisTopic models with different number of topics
cat(sprintf("Fit cisTopic models using runWarpLDAModels with %s topics. \n", k))
tmp.dir <- paste0(dirname(outfile), "/tmp/", project.name, ".", k, "topics")
dir.create(tmp.dir, showWarnings = F, recursive = T)

timing <- system.time({
  cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=k, nCores=nc, seed=seed,
                                     iterations = numiter, addModels=FALSE, tmp = tmp.dir)

})
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# SAVE RESULTS
# ------------
cat(sprintf("Saving results to %s \n", outfile))
saveRDS(list(k = k, fit = cisTopicObject),file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())


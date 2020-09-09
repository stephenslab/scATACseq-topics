#!/usr/bin/env Rscript

# Compile the fitted Poisson non-negative factorizations into a single
# .RData file.

library(tools)
library(stringr)
library(optparse)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser, c("--outdir","-o"),type="character",default="./")
parser <- add_option(parser, c("--dataset","-d"),type="character",default="dataset")
out    <- parse_args(parser)
out.dir  <- out$outdir
dataset  <- out$dataset
rm(parser,out)

# Combine results from all files of the form fit-*.rds in out.dir
# directory.

# List all the RDS files containing the model fits.
files <- Sys.glob(file.path(out.dir, paste0("fit-", dataset, "-*.rds")))
n     <- length(files)
cat(sprintf("Combined fitted results of %d files for %s dataset. \n", n, dataset))

# Set up two data structures: "fits", a list used to store all the
# results; and "dat", a data frame summarizing the model parameters
# and optimization settings used to produce these fits.
fits   <- vector("list",n)
labels <- files
labels <- str_remove(labels,paste0(out.dir,"/"))
labels <- str_remove(labels,".rds")
names(fits) <- labels
dat <- data.frame(label       = labels,
                  k           = 0,
                  method      = "em",
                  extrapolate = FALSE,
                  stringsAsFactors = FALSE)

# Load the results from the RDS files.
for (i in 1:n) {
  out                  <- readRDS(files[i])
  fits[[i]]            <- out$fit
  dat[i,"k"]           <- ncol(out$fit$F)
  dat[i,"method"]      <- out$method
  dat[i,"extrapolate"] <- out$control$extrapolate
}

# Reorder the results in "fits" and "dat".
dat  <- transform(dat,method = factor(method,c("em","ccd","scd")))
i    <- with(dat,order(k,extrapolate,method))
dat  <- dat[i,]
fits <- fits[i]
rownames(dat) <- NULL

# Convert the "k" column to a factor.
dat <- transform(dat,k = factor(k))

# Save the combined results to an .RData file.
out.file <- paste0(out.dir, "/compiled.fits.", dataset,".RData")
save(list = c("dat","fits"),
     file = out.file)
# resaveRdaFiles(out.file)

cat("Combined results saved to", out.file, "\n")


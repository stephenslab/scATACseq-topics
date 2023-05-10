# Here we convert the Poisson NMF to a binomial topic model for the
# Cusanovich et al (2018) data, with k = 13 topics. These
# were the steps taken to load R and allocate computing resources for
# this analysis:
#
#   sinteractive -p broadwl -c 8 --mem=16G --time=1:00:00
#   module load R/3.5.1
#

# Load a few packages.
library(tools)
library(Matrix)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data.
load("../data/Cusanovich_2018/processed_data/Cusanovich_2018.RData")

# Load the k = 10 Poisson NMF model fit.
fit_pois <-
  readRDS("../output/Cusanovich_2018/fit-Cusanovich2018-scd-ex-k=13.rds")$fit

# Convert the Poisson NMF model to a binomial topic model without any EM
# updates to refine the fit. (This step involves a simple rescaling of L
# and F, and should be very fast.)
fit_binom <- poisson2binom(counts,fit_pois,numem = 0)

# Save the results.
save(list = "fit_binom",file = "binom-fit-cusanovich2018-k=13.RData")
resaveRdaFiles("binom-fit-cusanovich2018-k=13.RData")

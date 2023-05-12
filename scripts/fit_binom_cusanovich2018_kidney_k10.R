# Here we convert the Poisson NMF to a binomial topic model for the
# Cusanovich et al (2018) data, kidney only, with k = 10 topics. These
# were the steps taken to load R and allocate computing resources for
# this analysis:
#
#   sinteractive -p broadwl -c 8 --mem=100G --time=8:00:00
#   module load R/3.5.1
#

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
fit_pois <- readRDS(
  file.path("../output/Cusanovich_2018/tissues",
            "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit

# Convert the Poisson NMF model to a binomial topic model without any EM
# updates to refine the fit. (This step involves a simple rescaling of L
# and F, and should be very fast.)
fit_binom <- poisson2binom(counts,fit_pois,numem = 0)

# Perform the conversion a second time, this time with several EM
# updates to refine the fit.
t0 <- proc.time()
fit_binom_em <- poisson2binom(counts,fit_pois,numem = 100)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results.
save(list = c("fit_binom","fit_binom_em"),
     file = "binom-fit-cusanovich2018-kidney-k=10.RData")
resaveRdaFiles("binom-fit-cusanovich2018-kidney-k=10.RData")

# Here we convert the Poisson NMF to a binomial topic model for the
# Buenrostro et al (2018) data, with k = 10 topics. These were the
# steps taken to load R and allocate computing resources for this
# analysis:
#
#   sinteractive -p broadwl -c 8 --mem=100G --time=1:00:00
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

# Convert the Poisson NMF model to a binomial topic model without any EM
# updates to refine the fit. (This step involves a simple rescaling of L
# and F, and should be very fast.)
fit_binom <- poisson2binom(counts,fit_pois,numem = 0)

# Perform the conversion a second time, this time with several EM
# updates to refine the fit.
t0 <- proc.time()
fit_binom_em <- poisson2binom(counts,fit_pois,numem = 20)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Compare the two binomial fits.
print(mean(abs(fit_binom$L - fit_binom_em$L) < 0.01))
print(mean(abs(fit_binom$L - fit_binom_em$L) < 0.05))
print(mean(abs(fit_binom$F - fit_binom_em$F) < 0.01))

# Save the results.
save(list = c("fit_binom","fit_binom_em"),
     file = "binom-fit-buenrostro2018-k=10.RData")
resaveRdaFiles("binom-fit-buenrostro2018-k=10.RData")

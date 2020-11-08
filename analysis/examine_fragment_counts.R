library(Matrix)
load("../data/Buenrostro_2018_counts.RData")

# Look at distribution of total counts across populations.
y <- rowSums(counts)/rowSums(counts > 0)
print(sort(tapply(y,factor(samples$label),mean)))

# Load the results.
dat <- readRDS("homer_motif_enrichment_results.rds")

# Compare the enriched motifs in the two MEP (erythroid) topics.
plot(dat$mlog10P[,"k5"],dat$mlog10P[,"k8"],pch = 20,
     xlab = "topic 5",ylab = "topic 8")

# Compare the enriched motifs in the two HSC/MPP topics.
plot(dat$mlog10P[,"k3"],dat$mlog10P[,"k9"],pch = 20,
     xlab = "topic 3",ylab = "topic 9")


# TO DO: Explain here what this script does, and how to use it.

# Load the results of the DE analysis using the k = 10 topic model.
source("code/motif_analysis.R")
selected_regions <- select_DA_regions(DA_res, method = "pval",
                                      thresh.pval = 0.01,
                                      out.dir = ".",
                                      save.bed = TRUE)

# For each topic, perform TF motif enrichment analysis using HOMER hypergeometric test.
cat("Performing motif enrichment analysis using HOMER.\n")
homer_res <- vector("list", ncol(DA_res$z))
names(homer_res) <- colnames(DA_res$z)
for(k in 1:ncol(DA_res$z)){
  homer_res[[k]] <- run_homer(selected_regions$filenames[k],
                              genome = genome,
                              homer.path = homerpath,
                              use.hypergeometric = TRUE,
                              out.dir=paste0(homer.dir, "/homer_result_topic_", k),
                              n.cores=nc)
}

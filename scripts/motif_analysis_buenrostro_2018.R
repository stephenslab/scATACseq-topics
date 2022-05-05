library(Matrix)
library(fastTopics)
source("code/motif_analysis.R")

# SETTINGS
# ---------
DAfile          <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/DAanalysis-Buenrostro2018-k=10/DA_regions_topics_noshrinkage_10000iters.rds'
genome          <- 'hg19'
homerpath       <- '/project2/xinhe/software/homer/bin/findMotifsGenome.pl'
nc              <- 8
out.dir         <- '/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/postfit_v2/motifanalysis-Buenrostro2018-k=10'

# DIFFERENTIAL ACCESSBILITY ANALYSIS
# ------------------------------------------

cat("Load precomputed differential accessbility statistics.\n")
DA_res <- readRDS(DAfile)

# Regions with NAs
rows_withNAs <- which(apply(DA_res$z, 1, anyNA))
cat(length(rows_withNAs), "regions with NAs in z-scores... \n")
head(DA_res$z[rows_withNAs,])


# SELECT REGIONS and PERFORM MOTIF ENRICHMENT ANALYSIS USING HOMER
# -----------------------------------------------------------------

# Select regions with DA p-value < 0.1
# -------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_pval_0.1_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "pval", thresh.pval = 0.1, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

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
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))


# Select regions with DA p-value < 0.05
# --------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_pval_0.05_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "pval", thresh.pval = 0.05, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

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
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))

# Select top 1% regions with largest logFC
# ------------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_top1percent_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "topPercent", top.percent = 0.01, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

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
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))

# Select top 2000 regions with largest logFC
# -------------------------------------------
homer.dir <- paste0(out.dir, "/HOMER/DA_top2000_regions")
cat("Select regions for motif enrichment analysis... \n")
cat(sprintf("%d regions in total. \n", nrow(DA_res$z)))

selected_regions <- select_DA_regions(DA_res, method = "topN", top.n = 2000, out.dir = homer.dir, save.bed = TRUE)
saveRDS(selected_regions, paste0(homer.dir, "/selected_regions.rds"))

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
saveRDS(homer_res, paste0(homer.dir, "/homer_knownResults.rds"))

sessionInfo()

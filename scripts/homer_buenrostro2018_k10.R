# TO DO: Explain here what this script does, and how to use it.
#
#   sinteractive -p broadwl ...
#   module load R/3.5.1
#
library(tools)
library(fastTopics)

# Load the results of the DE analysis using the k = 10 topic model.
load(file.path("../output/Buenrostro_2018/binarized/filtered_peaks",
               "de-buenrostro2018-k=10-noshrink.RData"))

# Get the feature positions.
feature_names <- rownames(de$postmean)
out           <- strsplit(feature_names,"_")
positions     <- data.frame(chr   = sapply(out,"[[",1),
                            start = sapply(out,"[[",2),
                            end   = sapply(out,"[[",3),
                            name  = feature_names,
                            stringsAsFactors = FALSE)

# Repeat for each topic.
k <- ncol(de$postmean)
for (i in 1:k) {

  # Create a BED file, "positions.bed", containing the regions with
  # p-value < 0.01.
  rows <- which(de$lpval[,i] > 2)
  write.table(positions[rows,],"positions.bed",sep = "\t",quote = FALSE,
              row.names = FALSE,col.names = FALSE)
}

source("../code/motif_analysis.R")
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

# Save the results.
# TO DO.

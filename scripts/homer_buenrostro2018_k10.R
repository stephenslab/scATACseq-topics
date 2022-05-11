# A short script for performing HOMER motif enrichment analysis using
# the differential accessibility results generated from the
# de_analysis_buenrostro2018_k10.R script. These were the steps taken
# to load R and allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 4 --mem=16G --time=12:00:00
#   module load R/3.5.1
#
# Also for this script to work the findMotifsGenome.pl Perl script
# needs to be findable (e.g., by adding the appropriate directory to
# the PATH environment variable).
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
topics <- colnames(de$postmean)
k      <- length(topics)
homer  <- vector("list",k)
names(homer) <- topics
for (i in 1:k) {

  # Create a BED file, "positions.bed", containing the regions with
  # p-value < 0.05.
  rows <- which(de$lpval[,i] > -log10(0.05))
  write.table(positions[rows,],"positions.bed",sep = "\t",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

  # Run the HOMER motif enrichment analysis.
  t0 <- proc.time()
  homer.command <- paste("findMotifsGenome.pl positions.bed hg19 homer",
                         "-len 8,10,12 -size 200 -mis 2 -S 25 -p 4 -h")
  system(homer.command)
  t1 <- proc.time()
  timing <- t1 - t0
  cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))
  homer[[i]] <- read.table("homer/knownResults.txt",sep = "\t",
                           comment.char = "",header = TRUE,
                           check.names = FALSE,stringsAsFactors = FALSE)
  system("rm -Rf positions.bed homer")
}

# Save the results.
saveRDS(homer,"homer-buenrostro2018-k=10-noshrink.rds")

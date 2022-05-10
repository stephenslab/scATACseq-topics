# TO DO: Explain here what this script does, and how to use it.
#
#   sinteractive -p broadwl -c 4 --mem=16G --time=24:00:00
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
  # p-value < 0.05.
  rows <- which(de$lpval[,i] > -log10(0.05))
  write.table(positions[rows,],"positions.bed",sep = "\t",quote = FALSE,
              row.names = FALSE,col.names = FALSE)

  # Run the HOMER motif enrichment analysis.
  homer.command <-
    paste("/scratch/midway2/pcarbo/homer/bin/findMotifsGenome.pl",
          "positions.bed hg19 homer -len 8,10,12 -size 200 -mis 2",
          "-S 25 -p 4 -h")
  system.out <- system(homer.command,ignore.stderr = TRUE,
                       ignore.stdout = TRUE,intern = TRUE)
  res <- read.table("homer/knownResults.txt",sep = "\t",comment.char = "",
                    header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
}

# Save the results.
# TO DO.

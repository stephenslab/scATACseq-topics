# TO DO: Explain here what this script is for, and how to use it.
library(Matrix)
library(tools)

# Load the previously prepared count data.
load(file.path("../data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))
ids <- rownames(counts)

# Load the gene activity scores downloaded from the Mouse sci-ATAC-seq
# Atlas website.
scores <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "activity_scores.quantitative.rds"))
scores <- scores["SLC12A1",ids]

# Save these data to an .RData file.
save(list = "scores",file = "slc12a1_data.RData")
resaveRdaFiles("slc12a1_data.RData")

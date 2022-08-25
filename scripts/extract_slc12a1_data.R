# Extract some data for the "explore_cicero" analysis.
library(Matrix)
library(tools)

# Load the previously prepared count data.
load(file.path("../data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))
ids <- rownames(counts)

# Get the gene activity scores downloaded from the Mouse sci-ATAC-seq
# Atlas website.
scores <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "activity_scores.quantitative.rds"))
scores <- scores["SLC12A1",ids]

# Get the peaks linked to gene Slc12a1 by Cicero, restricted to the
# kidney-related clusters (clusters 11, 18, 25 and 30).
cicero <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
cicero <- subset(cicero,
                 peak1.tss.gene_name == "Slc12a1" |
                 peak2.tss.gene_name == "Slc12a1")
cicero <- subset(cicero,is.element(cluster,c(11,18,25,30)))

# Save these data to an .RData file.
save(list = c("scores","cicero"),file = "slc12a1_data.RData")
resaveRdaFiles("slc12a1_data.RData")

# Extract some data for the "explore_cicero" analysis.
#
# The data files activity_scores.quantitative.rds and
# master_cicero_conns.rds files were downloaded from the Mouse
# sci-ATAC-seq Atlas website:
#
#   https://atlas.gs.washington.edu/mouse-atac/
#
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

# Get the peaks linked to some genes, restricted to the kidney-related
# clusters (clusters 11, 18, 25 and 30).
cicero <- readRDS(file.path("../data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
genes <- c("Slc24a5","Myef2","Ctxn2","Slc12a1","Dut","Fbn1","Cep152")
dat <- NULL
for (gene in genes)
  dat <- rbind(dat,
               subset(cicero,
                      peak1.tss.gene_name == gene |
                      peak2.tss.gene_name == gene))
dat <- subset(dat,is.element(cluster,c(11,18,25,30)))

# Get the peaks linked to gene Slc12a1 by Cicero, restricted to the
# kidney-related clusters (clusters 11, 18, 25 and 30).
cicero <- subset(dat,
                 peak1.tss.gene_name == "Slc12a1" |
                 peak2.tss.gene_name == "Slc12a1")
cicero_other <- dat

# Save these data to an .RData file.
save(list = c("scores","cicero","cicero_other"),
     file = "slc12a1_data.RData")
resaveRdaFiles("slc12a1_data.RData")

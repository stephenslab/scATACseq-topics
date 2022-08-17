# TO DO: Explain here what this script does, and how to use it.
#
# sinteractive -p mstephens --account=pi-mstephens -c 4 --mem=16G \
#   --time=24:00:00
# module load R/3.5.1
#
# Load a few packages.
library(tools)
library(Matrix)
library(cicero)
library(fastTopics)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the previously prepared count data, and reformat it as a
# "CellDataSet" (CDS) object.
load(file.path("../data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))
counts <- t(counts)
dat    <- summary(counts)
dat$i  <- rownames(counts)[dat$i]
dat$j  <- colnames(counts)[dat$j]
names(dat) <- c("peak","cell","count")
counts <- dat
rm(dat)
class(counts) <- "data.frame"
cds <- make_atac_cds(counts,binarize = TRUE)

# Compute low-dimensional embedding using t-SNE, and add the t-SNE
# co-ordinates to the CDS object.
cds <- detectGenes(cds)
cds <- reduceDimension(cds,max_components = 2,num_dim = 12,verbose = TRUE,
                       reduction_method = "tSNE",norm_method = "none")
tsne_coords <- t(reducedDimA(cds))
rownames(tsne_coords) <- rownames(pData(cds))
cicero_cds <- make_cicero_cds(cds,reduced_coordinates = tsne_coords)

stop()

t0 <- proc.time()
cons <- run_cicero(cicero_cds,sample_genome,sample_num = 2)
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

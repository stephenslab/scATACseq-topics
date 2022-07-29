library(Matrix)
library(fastTopics)
set.seed(1)
load("data/Cusanovich_2018/processed_data/Cusanovich_2018_metadata_only.RData")
fit <- readRDS("output/Cusanovich_2018/fit-Cusanovich2018-scd-ex-k=13.rds")$fit
x <- samples$cell_label
x[x == "Activated B cells"] <- "B cells"
x[x == "Immature B cells"] <- "B cells"
x[x == "Ex. neurons CPN"] <- "Neurons"
x[x == "Ex. neurons SCPN"] <- "Neurons"
x[x == "Ex. neurons CThPN"] <- "Neurons"
x[x == "SOM+ Interneurons"] <- "Neurons"
x[x == "Inhibitory neurons"] <- "Neurons"
x[x == "Regulatory T cells"] <- "T cells"
x[x == "Endothelial I cells"] <- "Endothelial cells"
x[x == "Endothelial II cells"] <- "Endothelial cells"
x[x == "Endothelial I (glomerular)"] <- "Endothelial cells"
x[x == "Proximal tubule S3"] <- "Proximal tubule"
x[x == "Type I pneumocytes"] <- "Pneumocytes"
x[x == "Type II pneumocytes"] <- "Pneumocytes"
x[x == "Alveolar macrophages"] <- "Macrophages"
# x[x == "Proximal tubule"] <- "cluster 18"
x[x == "Loop of henle"] <- "cluster 18"
x[x == "Distal convoluted tubule"] <- "cluster 18"
x[x == "Collecting duct"] <- "cluster 18"
x[x == "DCT/CD"] <- "cluster 18"
cell_types <-
  c("Cardiomyocytes",
    "Astrocytes",
    "Oligodendrocytes",
    "Hepatocytes",
    "Podocytes",
    "Endothelial cells",
    "Neurons",
    "Purkinje cells",
    "Cerebellar granule cells",
    "T cells",
    "NK cells",
    "B cells",
    "Dendritic cells",
    "Macrophages",
    "Microglia",
    "Monocytes",
    "Proximal tubule",
    "cluster 18",
    "Pneumocytes",
    "Enterocytes",
    "Erythroblasts",
    "Sperm",
    "Hematopoietic progenitors",
    "Collisions",
    "Unknown")
samples <- transform(samples,
                     tissue.replicate = factor(tissue.replicate),
                     cell_label = factor(x,cell_types))
topic_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
                  "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
                  "gray","darkblue","magenta")
pca_embed_method <- function (fit, ...)
    drop(pca_from_topics(fit,dims = 1))

# TO DO: Fix topic colours and order of tissues in first Structure plot.
p1 <- structure_plot(fit,grouping = samples$tissue.replicate,gap = 20,n = 4000,
                     perplexity = 30,colors = topic_colors,verbose = FALSE)
p2 <- structure_plot(fit,grouping = samples$cell_label,gap = 20,n = 4000,
                    perplexity = 30,colors = topic_colors,verbose = FALSE)
print(p2)

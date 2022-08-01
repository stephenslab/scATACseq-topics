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
                     cell_label = factor(x,cell_types))
p2 <- structure_plot(fit,grouping = samples$cell_label,gap = 20,n = 4000,
                    perplexity = 30,colors = topic_colors,verbose = FALSE)
print(p2)

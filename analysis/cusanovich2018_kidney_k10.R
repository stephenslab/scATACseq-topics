library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
load("data/Cusanovich_2018/processed_data/Cusanovich_2018_Kidney.RData")
fit <- readRDS(file.path("output/Cusanovich_2018/tissues",
                         "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit
fit <- poisson2multinom(fit)

# Plot using Cusanovich et al labeling.
set.seed(1)
x <- samples$cell_label
x[x == "Activated B cells" | x == "B cells" | x == "Alveolar macrophages" |
  x == "Monocytes" | x == "NK cells" | x == "T cells" | 
  x == "Regulatory T cells" | x == "Dendritic cells" |
  x == "Macrophages"] <- "Immune cells"
x[x == "Endothelial I (glomerular)" |
  x == "Endothelial I cells"] <- "Endothelial I cells"
x[x == "Collisions" | x == "Enterocytes" | x == "Hepatocytes" |
  x == "Sperm" | x == "Type II pneumocytes" | x == "Endothelial II cells" |
  x == "Hematopoietic progenitors" | x == "Unknown"] <- "Other or unknown"
x <- factor(x,c("Collecting duct","DCT/CD","Distal convoluted tubule",
                "Loop of henle","Proximal tubule","Proximal tubule S3",
                "Podocytes","Endothelial I cells","Immune cells",
                "Other or unknown"))
samples <- transform(samples,cell_label = x)
topic_colors <- c("peru","darkmagenta","magenta","gold","darkorange",
                  "red","lightgray","limegreen","dodgerblue","royalblue")
p1 <- structure_plot(fit,grouping = samples$cell_label,gap = 20,
                     colors = topic_colors,verbose = FALSE)
print(p1)

# Cluster the cells using the topic proportions only.
pca <- prcomp(fit$L)$x
pc1 <- pca[,1]
pc2 <- pca[,2]
x <- rep("other",nrow(pca))
x[pc1 < 0] <- "proximal tubule"
x[pc1 > 0 & pc2 > 0.2] <- "endothelial I"
rows <- which(x == "other")
fit2 <- select_loadings(fit,loadings = rows)
pca  <- prcomp(fit2$L)$x
y    <- rep("DCT + CD",nrow(pca))
pc1  <- pca[,1]
pc2  <- pca[,2]
y[pc1 < 0 & pc2 > -0.15] <- "loop of henle"
y[pc1 >= 0 & pc2 > -0.15] <- "podocytes + other"
x[rows] <- y
samples$cluster <- factor(x,c("DCT + CD","loop of henle","proximal tubule",
                              "podocytes + other","endothelial I"))

# Plot using our clustering.
set.seed(1)
p2 <- structure_plot(fit,grouping = samples$cluster,gap = 20,
                     colors = topic_colors,verbose = FALSE)
print(p2)

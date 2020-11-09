library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load the data.
load("../data/Buenrostro_2018_counts.RData")
samples <- transform(samples,label = factor(label))

# Load the Poisson NMF model fit.
fit <- readRDS("../output/fit-Buenrostro2018-scd-ex-k=11.rds")$fit

# Create PCA plots.
pca <- prcomp(poisson2multinom(fit)$L)$x
x   <- rep("A",nrow(pca))
pc1 <- pca[,1]
pc2 <- pca[,2]
pc3 <- pca[,3]
pc4 <- pca[,4]
pc5 <- pca[,5]
pc6 <- pca[,6]
pc8 <- pca[,8]
x[pc1 < 0.1]  <- "B"
x[pc2 < -0.2] <- "C"
x[pc3 > 0.35] <- "pDC"
x[pc4 > 0.55] <- "CLP"
x[pc5 < -0.25 & pc6 < -0.25] <- "mono"
x[pc8 > 0.35] <- "U"

# p1 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = samples$label) +
#   scale_fill_manual(values = topic_colors)
# p2 <- pca_plot(poisson2multinom(fit),pcs = 1:2,k = 5)
# p3 <- pca_hexbin_plot(poisson2multinom(fit),pcs = 1:2,bins = 30,
#                       breaks = c(0,1,5,10,20,Inf))
# p4 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = factor(x))
# plot_grid(p1,p2,p3,p4)

samples$cluster <- factor(x)
table(samples$label,samples$cluster)

# Project the cells---other than the "unknown" cluster---onto the top
# two PCs.
facs_colors <- c("skyblue",     # CLP
                 "gold",        # CMP
                 "orange",      # GMP
                 "forestgreen", # HSC
                 "turquoise",   # LMPP
                 "firebrick",   # MEP
                 "orangered",   # mono
                 "greenyellow", # MPP
                 "orchid",      # pDC
                 "gray")        # unknown
rows <- which(samples$cluster != "U")
fit2 <- select(poisson2multinom(fit),loadings = rows)
p1   <- pca_plot(fit2,pcs = 1:2,fill = samples[rows,"label"]) +
  scale_fill_manual(values = facs_colors) +
  labs(fill = "cell type")
p2   <- pca_plot(fit2,pcs = 3:4,fill = samples[rows,"label"]) +
  scale_fill_manual(values = facs_colors) +
  labs(fill = "cell type")
plot_grid(p1,p2)

# Create a Structure plot.
topic_colors <- c("greenyellow",
                  "yellow",      # CMP
                  "gold",        # CMP
                  "orchid",      # pDC
                  "orangered",   # mono
                  "firebrick",   # MEP
                  "orange",      # GMP
                  "skyblue",     # CLP
                  "gray",        # U
                  "forestgreen", # HSC + MPP
                  "turquoise")   # LMPP
topics <- c(2,9,8,4,5,3,11,1,7,6,10)
p3 <- structure_plot(poisson2multinom(fit),topics = topics,
                     colors = topic_colors[topics],
                     grouping = samples$cluster,gap = 20,
                     num_threads = 4,verbose = FALSE)
print(p3)

# Plot cluster B.
rows <- which(samples$cluster == "B")
fit2 <- select(poisson2multinom(fit),loadings = rows)
p4   <- pca_plot(fit2,fill = samples$label[rows,drop = FALSE]) +
  scale_fill_manual(values = facs_colors,drop = FALSE) +
  labs(fill = "cell type")

# Plot cluster C.
rows <- which(samples$cluster == "C")
fit2 <- select(poisson2multinom(fit),loadings = rows)
p5   <- pca_plot(fit2,fill = samples$label[rows,drop = FALSE]) +
  scale_fill_manual(values = facs_colors,drop = FALSE) +
  labs(fill = "cell type")


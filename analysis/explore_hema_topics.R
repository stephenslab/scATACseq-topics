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
topic_colors <- c("turquoise","darkorange","dodgerblue","gold","peru",
                  "greenyellow","firebrick","olivedrab","royalblue",
                  "forestgreen","gray")
pca <- prcomp(poisson2multinom(fit)$L)$x
x   <- rep("A",nrow(pca))
pc8 <- pca[,8]
x[pc8 > 0.35] <- "U"
p1 <- pca_plot(poisson2multinom(fit),pcs = c(7,8),fill = samples$label) +
  scale_fill_manual(values = topic_colors)
p2 <- pca_plot(poisson2multinom(fit),pcs = 7:8,k = 5)
p3 <- pca_hexbin_plot(poisson2multinom(fit),pcs = 7:8,bins = 30,
                      breaks = c(0,1,5,10,20,Inf))
p4 <- pca_plot(poisson2multinom(fit),pcs = 7:8,fill = factor(x))
plot_grid(p1,p2,p3,p4)

stop()

samples$cluster <- factor(x)

# Project the cells---other than the "unknown" cluster---onto the top
# two PCs.
facs_colors <- c("skyblue",     # CLP
                 "gold",        # CMP
                 "orange",      # GMP
                 "forestgreen", # HSC
                 "turquoise",   # LMPP
                 "firebrick",     # MEP
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
p6 <- structure_plot(poisson2multinom(fit),colors = topic_colors,
                     grouping = samples$label,gap = 20,
                     num_threads = 4)

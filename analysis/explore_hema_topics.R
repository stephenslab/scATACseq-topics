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
# fit <- readRDS("../output/fit-Buenrostro2018-binarized-scd-ex-k=10.rds")$fit
fit <- readRDS("../output/fit-Buenrostro2018-scd-ex-k=11.rds")$fit

# Create PCA plots.
topic_colors <- c("turquoise","darkorange","dodgerblue","gold","peru",
                  "greenyellow","firebrick","olivedrab","royalblue",
                  "forestgreen","gray")
p1 <- pca_plot(poisson2multinom(fit),pcs = 1:2,fill = samples$label) +
  scale_fill_manual(values = topic_colors)
p2 <- pca_hexbin_plot(poisson2multinom(fit),pcs = 1:2,bins = 30,
                      breaks = c(0,1,5,10,20,Inf))
plot_grid(p1,p2)

p3 <- pca_plot(poisson2multinom(fit),pcs = 3:4,fill = samples$label) +
  scale_fill_manual(values = topic_colors)
p4 <- pca_hexbin_plot(poisson2multinom(fit),pcs = 3:4,bins = 30,
                      breaks = c(0,1,5,10,20,Inf))
p5 <- pca_plot(poisson2multinom(fit),pcs = 3:4,k = 3)
plot_grid(p3,p4,p5,nrow = 1)

# Create a Structure plot.
p6 <- structure_plot(poisson2multinom(fit),colors = topic_colors,
                     grouping = samples$label,gap = 20,
                     num_threads = 4)

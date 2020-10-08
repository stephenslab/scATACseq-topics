# Short script to illustrate use of k-means to roughly define clusters
# based on the topic proportions, and use these clusters to create a
# Structure plot.
library(Matrix)
library(dplyr)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(1)

# Load data and Poisson NMF fit.
load("../data/Cusanovich_2018.RData")
fit <- readRDS("../output/fit-Cusanovich2018-scd-ex-k=13.rds")$fit

# Define clusters, roughly.
clusters <- factor(kmeans(poisson2multinom(fit)$L,centers = 14)$cluster)
print(sort(table(clusters),decreasing = TRUE))
      
# Create the Structure plot.
colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
            "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
            "gray")
n <- nrow(fit$L)
rows <- sample(n,4000)
p1 <- structure_plot(select(poisson2multinom(fit),loadings = rows),
                     grouping = clusters[rows],n = Inf,gap = 20,
                     perplexity = 50,topics = 1:13,colors = colors,
                     num_threads = 4,verbose = TRUE)
print(p1)

library(Matrix)
library(fastTopics)
set.seed(1)
load("data/Cusanovich_2018/processed_data/Cusanovich_2018_metadata_only.RData")
fit <- readRDS("output/Cusanovich_2018/fit-Cusanovich2018-scd-ex-k=11.rds")$fit
samples <- transform(samples,tissue.replicate = factor(tissue.replicate))
topic_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
                  "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928",
                  "gray","darkblue","magenta")
pca_embed_method <- function (fit, ...)
  drop(pca_from_topics(fit,dims = 1))
p <- structure_plot(fit,grouping = samples$tissue.replicate,gap = 20,n = 4000,
                    perplexity = 70,colors = topic_colors,verbose = FALSE,
                    embed_method = pca_embed_method)

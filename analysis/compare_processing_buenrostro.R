library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
# Data downloaded from original paper.
load("../data/Buenrostro_2018_binarized.RData")
# Load the topic model fit to the data downloaded from the original paper.
fit1 <- readRDS(file.path("../output/Buenrostro_2018",
                          "fit-Buenrostro2018-binarized-scd-ex-k=10.rds"))$fit
fit1 <- poisson2multinom(fit1)
# Load the topic model fit to data processed using Chen et al (2019) pipeline.
fit2 <- readRDS(file.path("../output/Buenrostro_2018_Chen2019pipeline",
                          "fit-Buenrostro2018-binarized-scd-ex-k=10.rds"))$fit
fit2 <- poisson2multinom(fit2)

labels <- factor(samples$label,c("mono","pDC","MEP","HSC","MPP","CLP",
                                 "LMPP","CMP","GMP","UNK"))
topic_colors <- c("darkorange","limegreen","magenta","gold","skyblue",
                  "darkblue","dodgerblue","darkmagenta","red","olivedrab")
set.seed(1)
p1 <- structure_plot(fit1,grouping = labels,colors = topic_colors,
                     gap = 20,perplexity = 50,verbose = TRUE) +
  ggtitle("data from original paper")
p2 <- structure_plot(fit2,grouping = labels,colors = topic_colors,
                     gap = 20,perplexity = 50,verbose = TRUE) +
  ggtitle("processed using Chen et al 2019 pipeline")
plot_grid(p1,p2,nrow = 2,ncol = 1)

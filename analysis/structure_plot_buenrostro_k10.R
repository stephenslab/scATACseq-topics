library(fastTopics)
library(Matrix)
set.seed(1)
load("Buenrostro_2018_binarized_counts.RData")
fit <- readRDS("../output/fit-Buenrostro2018-binarized-scd-k=10.rds")$fit
fit <- poisson2multinom(fit)
labels <- factor(samples$label,c("mono","pDC","MEP","HSC","MPP","CLP",
                                 "LMPP","CMP","GMP","UNK"))
topic_colors <- c("darkorange","limegreen","magenta","gold","skyblue",
                  "darkblue","dodgerblue","darkmagenta","red","olivedrab")
structure_plot(fit,grouping = labels,colors = topic_colors,
               gap = 20,perplexity = 50,verbose = FALSE)

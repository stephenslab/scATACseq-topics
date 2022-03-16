library(Matrix)
library(fastTopics)
library(ggplot2)
set.seed(1)
load("../data/Buenrostro_2018_binarized.RData")
labels <- factor(samples$label,c("mono","pDC","MEP","HSC","MPP","CLP",
                                 "LMPP","CMP","GMP","UNK"))
topic_colors <- c("darkorange","limegreen","magenta","gold","olivedrab",
                  "darkblue","dodgerblue","skyblue")
fit <- readRDS("fit-buenrostro-2018-k=8.rds")
fit <- poisson2multinom(fit)
p <- structure_plot(fit,grouping = samples$label,perplexity = 50,
                    colors = topic_colors,gap = 30,verbose = FALSE)
ggsave("structure_plot_buenrostro2018.png",p,height = 2,width = 8,
       dpi = 450,bg = "white")

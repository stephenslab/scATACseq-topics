library(Matrix)
library(fastTopics)
load("../data/Cusanovich_2018.RData")
samples$tissue <- factor(samples$tissue)
fit <- readRDS(file.path("../output/Cusanovich_2018",
                         "fit-Cusanovich2018-scd-ex-k=10.rds"))$fit
fit <- poisson2multinom(fit)
topic_colors <- c("darkorange","limegreen","magenta","gold","skyblue",
                  "darkblue","dodgerblue","darkmagenta","red","olivedrab")
set.seed(1)
p <- structure_plot(fit,n = 2000,grouping = samples$tissue,
                    colors = topic_colors,gap = 30,perplexity = 100,
                    verbose = FALSE)
print(p)


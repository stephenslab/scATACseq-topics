library(Matrix)
library(fastTopics)
library(ggplot2)
set.seed(1)

outdir <- "/project2/mstephens/kevinluo/scATACseq-topics/output/Buenrostro_2018/binarized/filtered_peaks"

load(file.path(outdir, "Buenrostro_2018_binarized_filtered.RData"))
labels <- factor(samples$label,c("mono","pDC","MEP","HSC","MPP","CLP",
                                 "LMPP","CMP","GMP","UNK"))
topic_colors <- c("darkorange","limegreen","magenta","gold","olivedrab",
                  "darkblue","dodgerblue","skyblue")
fit <- readRDS(file.path(outdir, "fit-Buenrostro2018-binarized-k=8"))

# Load the K = 10 topic model fit to filtered data
fit <- readRDS(file.path(outdir, "/fit-Buenrostro2018-binarized-k=8"))$fit

fit <- poisson2multinom(fit)
p <- structure_plot(fit,grouping = samples$label,perplexity = 50,
                    colors = topic_colors,gap = 30,verbose = FALSE)
p
ggsave(file.path(outdir, "structure_plot_buenrostro2018_filtered.png"),p,height = 2,width = 8,
       dpi = 450,bg = "white")

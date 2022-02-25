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
p <- structure_plot(fit,n = 4000,grouping = samples$tissue,
                    colors = topic_colors,gap = 30,perplexity = 100,
                    verbose = FALSE)
print(p)

# Look closely at topics explaining chromatin accessibility in brain
# tissues (cerebellum, prefrontal cortex, whole brain).
brain_tissues <- c("Cerebellum","PreFrontalCortex","WholeBrain")
i <- which(is.element(samples$tissue,brain_tissues))
j <- which(!is.element(samples$tissue,brain_tissues))
plot(colMeans(fit$L[j,]),colMeans(fit$L[i,]),pch = 20,col = topic_colors,
     xlim = c(0,0.5),ylim = c(0,0.5),xlab = "non-brain cells",
     ylab = "brain cells")
abline(a = 0,b = 1,col = "gainsboro",lty = "dotted")

# Look closely at the topics explaining chromatin accessibility in
# kidney cells.
i <- which(is.element(samples$tissue,"Kidney"))
j <- which(!is.element(samples$tissue,"Kidney"))
plot(colMeans(fit$L[j,]),colMeans(fit$L[i,]),pch = 20,col = topic_colors,
     xlim = c(0,0.6),ylim = c(0,0.6),xlab = "non-kidney cells",
     ylab = "kidney cells")
abline(a = 0,b = 1,col = "gainsboro",lty = "dotted")

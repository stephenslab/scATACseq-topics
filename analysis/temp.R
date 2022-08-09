library(Matrix)
library(fastTopics)
load(file.path("output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))
k <- "k8"
plot(de$postmean[,k],pmin(40,de$lpval[,k]),pch = 20,
     xlab = "LFC",ylab = "-log10pval")

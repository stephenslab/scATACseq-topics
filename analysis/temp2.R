library(Matrix)
k <- "k8"
load(file.path("data/Cusanovich_2018/processed_data",
               "Cusanovich_2018_Kidney.RData"))
cicero <- readRDS(file.path("data/Cusanovich_2018/processed_data",
                            "master_cicero_conns.rds"))
cicero_slc12a1 <- subset(cicero,
                         peak1.tss.gene_name == "Slc12a1" |
                         peak2.tss.gene_name == "Slc12a1")
peaks <- unique(c(as.character(cicero_slc12a1$Peak1),
                  as.character(cicero_slc12a1$Peak12)))
load(file.path("output/Cusanovich_2018/tissues",
               "de-cusanovich2018-kidney-k=10-noshrink.RData"))
plot(de$postmean[peaks,k],pmin(40,de$lpval[peaks,k]),
     pch = 20,xlab = "LFC",ylab = "-log10pval",main = k)

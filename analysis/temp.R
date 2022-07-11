library(ggplot2)
library(cowplot)

# Draft analysis of motifs in topics k = 1 through 10.
homer <- readRDS(file.path("../output/Buenrostro_2018",
                           "binarized/filtered_peaks",
                           "homer-buenrostro2018-k=10-noshrink.rds"))

# Compile the motif enrichment p-values into a single table, after
# removing duplicate results.
motifs <- sort(unique(homer$k1[,"Motif Name"]))
n      <- length(motifs)
pvals  <- matrix(0,n,k)
rownames(pvals) <- motifs
colnames(pvals) <- paste0("k",1:k)
for (i in 1:k) {
  dat  <- homer[[i]]
  rows <- which(!duplicated(dat[,"Motif Name"]))
  dat  <- dat[rows,]
  rownames(dat) <- dat[,"Motif Name"]
  pvals[,i] <- dat[motifs,"P-value"]
}
pvals <- as.data.frame(pvals)

# For each topic, pick the top 10 motifs by p-value.
motifs <- NULL
for (i in 1:k) {
  rows1  <- which(pvals[,i] <= 0)
  rows2  <- order(pvals[,i])
  motifs <- c(motifs,rownames(pvals)[c(rows1,rows2[1:10])])
}
motifs <- sort(unique(motifs))
pvals <- pvals[motifs,]

# Plot the p-values in a tile plot.
# Colors from colorbrewer2.org.
colors <- c("#d73027","#f46d43","#fdae61","#fee090",
            "#e0f3f8","#abd9e9","#74add1",",#4575b4")
pdat <- NULL
for (i in 1:k) {
  pdat <- rbind(pdat,
                data.frame(motif = motifs,
                           topic = paste0("k",i),
                           lpval = log10(pvals[,i] + 1e-256)))
}
pdat <- transform(pdat,
                  lpval = cut(lpval,seq(-257,0,length.out = 8)),
                  motif = factor(motif,rev(motifs)))
p <- ggplot(pdat,aes(x = topic,y = motif,fill = lpval)) +
  geom_tile(color = "white",size = 0.5) +
  scale_fill_manual(values = colors) +
  theme_cowplot(font_size = 8)

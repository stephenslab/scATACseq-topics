library(fastTopics)
library(pathways)
library(ggplot2)
library(cowplot)

# Load the base-pair positions of the genes for the mm9 assembly.
load("data/mm9_seq_gene.RData")

# Load the data on gene Slc12a1.
load("data/Cusanovich_2018/processed_data/slc12a1_data.RData")

# Load the K = 10 topic model fit.
fit <- readRDS(file.path("output/Cusanovich_2018/tissues",
                         "fit-Cusanovich2018-Kidney-scd-ex-k=10.rds"))$fit
fit <- poisson2multinom(fit)

# Load the results of the DE analysis using the k = 10 topic model,
# without the adaptive shrinkage step.
load(file.path("output/Cusanovich_2018/tissues",
	       "de-cusanovich2018-kidney-k=10-noshrink.RData"))

# Get the base-pair positions of the peaks.
feature_names <- rownames(de$postmean)
out           <- strsplit(feature_names,"_")
positions     <- data.frame(chr   = sapply(out,"[[",1),
                            start = sapply(out,"[[",2),
                            end   = sapply(out,"[[",3),
                            name  = feature_names,
                            stringsAsFactors = FALSE)
positions <- transform(positions,
                       start = as.numeric(start),
                       end   = as.numeric(end))

# This plot shows that there is a strong relationship between the
# "gene activity scores" and the topic proportions for topic 8.
pdat <- data.frame(loading = fit$L[,"k8"],score = log10(1 + scores))
b <- coef(lm(score ~ loading,data = pdat))
p1 <- ggplot(pdat,aes(x = loading,y = score)) +
  geom_point() +
  geom_abline(intercept = b["(Intercept)"],slope = b["loading"],
              color = "dodgerblue",size = 1,linetype = "dashed") +
  labs(x = "topic proportion",y = "gene activity score") +
  theme_cowplot()

# Here let's examine the de_analysis results for the peaks near gene
# Slc12a1.
d <- 4e5
load("data/mm9_seq_gene.RData")
seq_gene <- subset(seq_gene,feature_name == "Slc12a1")
rows <- which(positions$chr == paste0("chr",seq_gene$chromosome) &
              positions$start > seq_gene$chr_start - d &
              positions$end < seq_gene$chr_stop + d)
peaks <- c(as.character(cicero$Peak1),as.character(cicero$Peak2))

x <- factor(is.element(positions[rows,"name"],peaks))
levels(x) <- c("dodgerblue","darkorange")
x <- as.character(x)
plot(positions[rows,"start"],de$lpval[rows,"k8"],pch = 20,col = x)
plot(positions[rows,"start"],de$postmean[rows,"k8"],pch = 20)

library(fastTopics)
library(pathways)
library(ashr)
library(ggplot2)
library(cowplot)

# Load the base-pair positions of the genes for the mm9 assembly.
load("data/mm9_seq_gene.RData")

# Load the data on gene Slc12a1.
load("data/Cusanovich_2018/processed_data/slc12a1_data.RData")
cicero <- transform(cicero,
                    Peak1 = as.character(Peak1),
                    Peak2 = as.character(Peak2))

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
# Slc12a1; specifically, within 400 kb of the transcribed region of
# gene Slc12a1.
d <- 4e5
seq_gene <- subset(seq_gene,feature_name == "Slc12a1")
rows <- which(with(positions,chr == paste0("chr",seq_gene$chromosome) &
                   start > seq_gene$chr_start - d &
                   end < seq_gene$chr_stop + d))
peaks <- c(cicero$Peak1,cicero$Peak2)
pdat <- data.frame(start      = positions[rows,"start"]/1e5,
                   postmean   = de$postmean[rows,"k8"],
                   z          = de$z[rows,"k8"],
                   lpval      = de$lpval[rows,"k8"],
                   cicero_hit = is.element(positions[rows,"name"],peaks))
p1 <- ggplot(pdat,aes(x = start,y = lpval,color = cicero_hit)) +
  geom_point() +
  scale_color_manual(values = c("darkblue","magenta")) +
  theme_cowplot()
p2 <- ggplot(pdat,aes(x = start,y = postmean,color = cicero_hit)) +
  geom_point() +
  geom_hline(yintercept = 0,color = "black",linetype = "dashed") +
  scale_color_manual(values = c("darkblue","magenta")) +
  theme_cowplot()
plot_grid(p1,p2,nrow = 2,ncol = 1)

# Run ashr on these results.
pdat <- transform(pdat,se = postmean/z)
fit <- ash(pdat$postmean,pdat$se,mixcompdist = "normal",method = "shrink")
pdat$ash <- fit$result$PosteriorMean
ggplot(pdat,aes(x = postmean,y = ash,color = se)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "black",linetype = "dashed") +
  theme_cowplot()

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

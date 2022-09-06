# Slc24a5, Myef2, Ctxn2, Slc12a1, Dut, Fbn1
gene <- "Dut"
class(cicero) <- "data.frame"
class(cicero_other) <- "data.frame"
dat <- subset(cicero_other,
              peak1.tss.gene_name == gene |
              peak2.tss.gene_name == gene)
peaks <- unique(c(as.character(dat$Peak1),
                  as.character(dat$Peak2)))
print(length(peaks))
pdat <- data.frame(start      = positions[rows,"start"]/1e5,
                   postmean   = de$postmean[rows,k],
                   z          = de$z[rows,k],
                   lpval      = de$lpval[rows,k],
                   cicero_hit = is.element(positions[rows,"name"],peaks))
p1 <- ggplot(pdat,aes(x = start,y = lpval,color = cicero_hit)) +
  geom_point() +
  scale_color_manual(values = c("royalblue","magenta")) +
  labs(x = "position on chromosome 2 (kb)",
       y = "-log10 p-value") +
  theme_cowplot()
p2 <- ggplot(pdat,aes(x = start,y = postmean,color = cicero_hit)) +
  geom_point() +
  geom_hline(yintercept = 0,color = "black",linetype = "dashed") +
  scale_color_manual(values = c("royalblue","magenta")) +
  labs(x = "position on chromosome 2 (kb)",
       y = "LFC estimate") +
  theme_cowplot()
print(plot_grid(p1,p2,nrow = 2,ncol = 1))

dat <- subset(cicero_other,
              Peak1 == "chr2_124977295_124978483" |
              Peak2 == "chr2_124977295_124978483")
print(unique(c(as.character(dat$peak1.tss.gene_name),
               as.character(dat$peak2.tss.gene_name))))

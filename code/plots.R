# For each topic, plot the number of samples exceeding specified topic
# proportion, as specified by "probs". The topics are ordered in the
# bar chart from most abundant to least abundant.
create_abundance_plot <- function (fit, probs = c(0.1,0.25,0.5),
                                   clrs=c("tomato","dodgerblue","darkblue")) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  m      <- length(probs)
  k      <- ncol(fit$L)
  L      <- fit$L
  out    <- apply(L,2,function (x) rev(cumsum(rev(table(cut(x,c(probs,1)))))))
  topics <- colnames(L)
  topics <- topics[order(out[1,],decreasing = TRUE)]
  dat    <- data.frame(topic = factor(rep(colnames(L),each = m),topics),
                       prob  = factor(rep(probs,times = k)),
                       n     = as.vector(out))
  return(ggplot(dat,aes_string(x = "topic",y = "n",fill = "prob")) +
           geom_col(color = "white",position = "dodge",width = 0.8) +
           scale_fill_manual(values = clrs) +
           labs(x = "topic",y = "samples") +
           theme_cowplot(font_size = 10))
}

# Create a basic scatterplot showing the topic proportions projected
# onto two principal components (PCs), and the colour of the points is
# varied according to a factor ("labels").
pca_plot_with_labels <-
  function (fit, pcs = c("PC1","PC2"), labels,
            colors = c("darkorange","darkblue","dodgerblue","magenta",
                       "gold")) {
    if (inherits(fit,"poisson_nmf_fit"))
      fit <- poisson2multinom(fit)
    out.pca <- prcomp(fit$L)
    dat     <- cbind(out.pca$x,data.frame(label = factor(labels)))
    return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2],fill = "label")) +
             geom_point(shape = 21,color = "white",size = 1.25) +
             scale_fill_manual(values = colors) +
             theme_cowplot(font_size = 10))
  }

# Create a "hex plot" showing the density of the data points
# (specifically, the topic proportions) as they are projected onto two
# principal ccomponents (PCs).
pca_hex_plot <-
  function (fit, pcs = c("PC1","PC2"), n = 40, bins = c(0,1,10,100,1000,Inf),
            colors = c("gainsboro","lightskyblue","gold","orange","magenta")) {
    if (inherits(fit,"poisson_nmf_fit"))
      fit <- poisson2multinom(fit)
    pdat <- as.data.frame(prcomp(fit$L)$x)
    return(ggplot(pdat,aes_string(x = pcs[1],y = pcs[2])) +
             stat_bin_hex(mapping = aes_q(fill = quote(cut(..count..,bins))),
                          bins = n) +
             scale_fill_manual(values = colors) +
             labs(fill = "count") +
             theme_cowplot(font_size = 10))
  }

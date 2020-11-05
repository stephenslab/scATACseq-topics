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
# onto two principal components (PCs).
basic_pca_plot <- function (fit, pcs = 1:2) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(prcomp(fit$L)$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
         geom_point(shape = 21,color = "white",fill = "black",size = 1.25) +
         theme_cowplot(font_size = 10))
}

# Create a basic scatterplot showing the topic proportions projected
# onto two principal components (PCs), and the colour of the points is
# varied according to a factor ("labels").
labeled_pca_plot <-
  function (fit, pcs = 1:2, labels, font_size = 9,
            colors = rep(c("firebrick","dodgerblue","forestgreen",
                           "darkmagenta","darkorange","gold","darkblue",
                           "peru","greenyellow","olivedrab"),times = 4),
            shapes = rep(c(21:24),each = 10),
            legend_label = "label") {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(prcomp(fit$L)$x)
  if (is.numeric(pcs))
    pcs <- names(dat)[pcs]
  dat <- cbind(data.frame(label = factor(labels)),dat)
  return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2],fill = "label",
                               shape = "label")) +
         geom_point(color = "white",size = 1.2,na.rm = TRUE) +
         scale_fill_manual(values = colors) +
         scale_shape_manual(values = shapes) +
         labs(fill = legend_label,shape = legend_label) +
         theme_cowplot(font_size = font_size))
}

# # Create a "hexbin plot" showing the density of the data points
# # (specifically, the topic proportions) as they are projected onto two
# # principal components (PCs).
# pca_hexbin_plot <-
#   function (fit, pcs = 1:2, n = 40, bins = c(0,1,10,100,1000,Inf),
#             colors = c("gainsboro","lightskyblue","gold","orange","magenta")) {
#   if (inherits(fit,"poisson_nmf_fit"))
#     fit <- poisson2multinom(fit)
#   dat <- as.data.frame(prcomp(fit$L)$x)
#   if (is.numeric(pcs))
#     pcs <- names(dat)[pcs]
#   return(ggplot(dat,aes_string(x = pcs[1],y = pcs[2])) +
#          stat_bin_hex(mapping = aes_q(fill = quote(cut(..count..,bins))),
#                       bins = n) +
#          scale_fill_manual(values = colors) +
#          labs(fill = "count") +
#          theme_cowplot(font_size = 10))
# }

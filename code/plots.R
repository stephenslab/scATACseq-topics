# Settings for GSEA plots
gsea_plot_shapes <- c(23,24,21)
gsea_plot_colors <- c("gold",         # biocyc
                      "yellowgreen",  # C1
                      "dodgerblue",   # C2
                      "olivedrab",    # C3
                      "firebrick",    # C4
                      "darkorange",   # C5
                      "magenta",      # C6
                      "darkblue",     # C7
                      "gold",         # H
                      "tomato",       # humancyc
                      "olivedrab",    # inoh
                      "darkblue",     # kegg
                      "royalblue",    # netpath
                      "darkorange",   # panther
                      "yellowgreen",  # pathbank
                      "magenta",      # pid
                      "dodgerblue")   # reactome

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

# Create a basic density plot for the region probability for each topic.
basic_density_plot <- function (fit, k = 1) {
  if (inherits(fit,"poisson_nmf_fit"))
    fit <- poisson2multinom(fit)
  dat <- as.data.frame(fit$F)
  if (is.numeric(k))
    name_topic <- names(dat)[k]
  return(ggplot(dat, aes_string(x=name_topic)) +
           geom_density(color="black", fill="white") +
           labs(x = "Probability", y = "Density", title = paste("Topic",k)) +
           theme_cowplot(font_size = 10))
}


#' @rdname genescore_volcano_plot
#'
#' @title Volcano plot of gene level statistics
#'
#' @description Create one or more "volcano" plots to visualize the
#'   results of gene level statistics computed from region level
#'   statistics using a topic model.
#'
#' @details Only points above a specified z-score quantile are labeled.
#' Note that genes with a weighted mean accessbility of zero are not shown.
#'
#' To better accommodate situations in which some z-scores
#' are much larger than all the others, the z-scores
#' are plotted on the square-root scale. To change
#' this, as well as other aspects, of the plot, replace
#' \code{volcano_plot_ggplot_call} with your own function; see input
#' argument \dQuote{ggplot_call}.
#'
#' The \dQuote{ggrepel} package is used to arrange the labels in a
#' visually attractive manner.
#'
#' Use interactive volcano plot is created using the \dQuote{plotly}
#' package. The "hover text" shows the label (see input argument
#' \dQuote{labels}) and detailed log-fold change statistics.
#'
#' @param genescore_res A list of gene level statistics computed from
#'   a distance weighted model of region level statistics
#'   obtained using \code{\link{diff_count_analysis}}.
#'
#' @param k The topic, or topics, selected by number or name. When not
#'   specified, all topics are plotted.
#'
#' @param title Figure title.
#'
#' @param labels Character vector specifying how the points in the
#'   volcano plot are labeled. This should be a character vector with
#'   one entry per log-fold change statistic (row of
#'   \code{genescore_res$beta}). When not specified, the row names
#'   of \code{genescore_res$beta} are used, if available. Labels are
#'   added to the plot using \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param betamax Truncate the log-fold change statistics
#'   by this amount. Any statistics greater in magnitude
#'   than \code{betamax} are set to \code{betamax} or \code{-betamax}.
#'
#' @param label_above_quantile Only z-scores above this quantile are
#'   labeled in the volcano plot. \code{\link[ggrepel]{geom_text_repel}}
#'   will attempt to label  all points when \code{label_above_quantile = 0}.
#'   When \code{label_above_quantile = Inf}, no points are labeled.
#'
#' @param subsample_below_quantile A number between 0 and 1. If
#'   greater than zero, log-fold change statistics with z-scores
#'   below this quantile will be subsampled according to
#'   \code{subsample_rate}. This is useful for large data sets to to
#'   reduce the number of points plotted.
#'
#' @param subsample_rate A number between 0 and 1 giving the
#'   proportion of log-fold change statistics with "small" z-scores that
#'   are included in the plot, uniformly at random. This is only used if
#'   \code{subsample_below_quantile} is greater than zero.
#'
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param ggplot_call The function used to create the plot. Replace
#'   \code{volcano_plot_ggplot_call} with your own function to customize
#'   the appearance of the plot.
#'
#' @param plot_grid_call When multiple topics are selected, this is
#'   the function used to arrange the plots into a grid using
#'   \code{\link[cowplot]{plot_grid}}. It should be a function accepting
#'   a single argument, \code{plots}, a list of \code{ggplot} objects.
#'
#' @return A \code{ggplot} object or a \code{plotly} object.
#'
#' @importFrom cowplot plot_grid
#'
#' @export
#'
genescore_volcano_plot <-
  function (genescore_res, k, labels,
            betamax = 10, label_above_quantile = 0.99,
            subsample_below_quantile = 0, subsample_rate = 0.1,
            max.overlaps = Inf, ggplot_call = genescore_volcano_plot_ggplot_call,
            plot_grid_call = function (plots) do.call(plot_grid,plots)) {

    # Check and process input arguments.
    beta <- genescore_res$beta

    if (missing(k))
      k <- seq(1,ncol(beta))

    if (missing(labels)) {
      if (!is.null(rownames(beta)))
        labels <- rownames(beta)
      else
        labels <- as.character(seq(1,nrow(beta)))
    }

    if (!(is.character(labels) & length(labels) == nrow(beta)))
      stop("Input argument \"labels\", when specified, should be a character ",
           "vector with one entry per log-fold change statistic (column of ",
           "the counts matrix)")

    if (length(k) == 1) {
      y.label <- "|z-score|"
      dat <- compile_genescore_volcano_plot_data(genescore_res,k,labels,
                                                 betamax, label_above_quantile,
                                                 subsample_below_quantile,subsample_rate)
      return(ggplot_call(dat,y.label,k,max.overlaps))
    } else {
      # Create a volcano plot for each selected topic, and combine them
      # using plot_grid. This is done by recursively calling genescore_volcano_plot.
      m     <- length(k)
      plots <- vector("list",m)
      names(plots) <- k
      for (i in 1:m)
        plots[[i]] <- genescore_volcano_plot(genescore_res,k[i],labels,betamax,
                                             label_above_quantile,subsample_below_quantile,
                                             subsample_rate,max.overlaps,ggplot_call,NULL)
      return(plot_grid_call(plots))

    }
  }

#' @rdname genescore_volcano_plot
#'
#' @param file Save the interactive volcano plot to this HTML
#'   file using \code{\link[htmlwidgets]{saveWidget}}.
#'
#' @param width Width of the plot in pixels. Passed as argument
#'   \dQuote{width} to \code{\link[plotly]{plot_ly}}.
#'
#' @param height Height of the plot in pixels. Passed as argument
#'   \dQuote{height} to \code{\link[plotly]{plot_ly}}.
#'
#' @param title The text used for the plot title.
#'
#' @param plot_ly_call The function used to create the plot. Replace
#'   \code{genescore_volcano_plot_ly_call} with your own function to customize
#'   the appearance of the interactive plot.
#'
#' @importFrom htmlwidgets saveWidget
#'
#' @export
#'
genescore_volcano_plotly <- function (genescore_res, k, file, labels,
                                      betamax = 10,
                                      subsample_below_quantile = 0,
                                      subsample_rate = 0.1, width = 600, height = 500,
                                      title = paste("topic",k),
                                      plot_ly_call = genescore_volcano_plot_ly_call) {

  # Check and process input arguments.
  beta <- genescore_res$beta

  if (missing(labels)) {
    if (!is.null(rownames(beta)))
      labels <- rownames(beta)
    else
      labels <- as.character(seq(1,nrow(beta)))
  }

  if (!(is.character(labels) & length(labels) == nrow(beta)))
    stop("Input argument \"labels\", when specified, should be a character ",
         "vector with one entry per log-fold change statistic (column of ",
         "the counts matrix)")

  # Compile the plotting data.
  y.label <- "|z-score|"
  dat <- compile_genescore_volcano_plot_data(genescore_res,k,labels,
                                             betamax, label_above_quantile = 0,
                                             subsample_below_quantile,subsample_rate)

  # Create the interactive volcano plot using plotly.
  p <- plot_ly_call(dat,y.label,title,width,height)
  if (!missing(file))
    saveWidget(p,file,selfcontained = TRUE,title = title)
  return(p)
}

# This is used by genescore_volcano_plot to compile the gene score results
# passed to ggplot.
#
#' @importFrom stats quantile
compile_genescore_volcano_plot_data <- function (genescore_res, k, labels,
                                                 betamax = 10, label_above_quantile = 0.99,
                                                 subsample_below_quantile = 0,
                                                 subsample_rate = 0.1) {

  dat <- with(genescore_res,
              data.frame(label = labels,
                         mean  = colmeans,
                         beta  = beta[,k],
                         z     = Z[,k],
                         y     = 0,
                         stringsAsFactors = FALSE))

  dat$y <- abs(genescore_res$Z[,k])
  rows     <- which(dat$mean > 0)
  dat      <- dat[rows,]
  dat <- transform(dat,beta = sign(beta) * pmin(betamax,abs(beta)))
  if (is.infinite(label_above_quantile))
    y0 <- Inf
  else
    y0 <- quantile(dat$y,label_above_quantile)
  dat$label[dat$y < y0] <- ""
  if (subsample_below_quantile > 0) {
    y0    <- quantile(dat$y,subsample_below_quantile)
    rows1 <- which(dat$y >= y0)
    rows2 <- which(dat$y < y0)
    rows2 <- sample(rows2,ceiling(subsample_rate * length(rows2)))
    rows  <- sort(c(rows1,rows2))
    message(sprintf("%d out of %d data points will be included in plot",
                    length(rows),nrow(dat)))
    dat   <- dat[rows,]
  }
  return(dat)
}

#' @rdname genescore_volcano_plot
#'
#' @param dat A data frame passed as input to \code{\link[ggplot2]{ggplot}}.
#'
#' @param y.label Label to use in the plot for \dQuote{dat$y}.
#'
#' @param topic.label The name or number of the topic being plotted.
#'   Only used to determine the plot title.
#'
#' @param max.overlaps Argument passed to
#'   \code{\link[ggrepel]{geom_text_repel}}.
#'
#' @param font.size Font size used in plot.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
genescore_volcano_plot_ggplot_call <- function (dat, y.label, topic.label,
                                                max.overlaps = Inf, font.size = 9) {
  ggplot(dat,aes_string(x = "beta",y = "y",fill = "mean",label = "label")) +
    geom_point(color = "white",stroke = 0.3,shape = 21,na.rm = TRUE) +
    scale_y_continuous(trans = "sqrt",
                       breaks = c(0,1,2,5,10,20,50,100,200,500,1e3,2e3,5e3,1e4,2e4,5e4)) +
    scale_fill_gradient2(low = "deepskyblue",mid = "gold",high = "orangered",
                         na.value = "gainsboro",
                         midpoint = mean(range(dat$mean))) +
    geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                    segment.color = "black",segment.size = 0.25,
                    max.overlaps = max.overlaps,na.rm = TRUE) +
    labs(x = "log-fold change",y = y.label,fill = "mean",
         title = paste("topic",topic.label)) +
    theme_cowplot(font.size)
}

#' @rdname genescore_volcano_plot
#'
#' @importFrom plotly plot_ly
#' @importFrom plotly hide_colorbar
#' @importFrom plotly layout
#'
#' @export
#'
genescore_volcano_plot_ly_call <- function (dat, y.label, title, width, height) {
  p <- plot_ly(data = dat,x = ~beta,y = ~sqrt(y),color = ~mean,
               colors = c("deepskyblue","gold","orangered"),
               text = ~sprintf(paste0("%s\nmean: %0.3f\n\u03b2: %+0.3f\n",
                                      "z: %+0.3f"),
                               label,mean,beta,z),
               type = "scatter",mode = "markers",hoverinfo = "text",
               width = width,height = height,
               marker = list(line = list(color = "white",width = 1),size = 7.5))
  p <- hide_colorbar(p)
  p <- layout(p,xaxis = list(title = "log-fold change (\u03b2)",
                             zeroline = FALSE,showgrid = FALSE),
              yaxis = list(title = paste("sqrt",y.label),
                           zeroline = FALSE,showgrid = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              showlegend = FALSE,title = title)
  return(p)
}

# Compile the data frame used for create_gsea_plot and create_gsea_plotly.
compile_data_for_gsea_plot <- function (gene_set_info, gsea_res, k,
                                        max_name_len = 44) {

  # Compute the "signed -log10(p-value)".
  P <- with(gsea_res,-sign(ES) * log10(pval))

  # Compute the "most extreme" signed -log10(p-value) among other topics.
  n  <- nrow(P)
  p0 <- rep(0,n)
  for (i in 1:n) {
    if (is.na(P[i,k]))
      p0[i] <- NA
    else if (sign(gsea_res$ES[i,k]) > 0)
      p0[i] <- max(P[i,-k],na.rm = TRUE)
    else
      p0[i] <- min(P[i,-k],na.rm = TRUE)
  }

  # Remove the non-breaking space character in some name of gene_set_info
  gene_set_info$name <- gsub("\xa0"," ",gene_set_info$name)

  # Compile the data for plotting.
  dat <- data.frame(p1         = P[,k],
                    p0         = p0,
                    name       = substr(gene_set_info$name,1,max_name_len),
                    id         = gene_set_info$id,
                    database   = gene_set_info$database,
                    collection = with(gene_set_info,
                                      ifelse(database == "MSigDB",
                                             as.character(category_code),
                                             as.character(data_source))),
                    stringsAsFactors = FALSE)

  # Subsample the gene sets with -log10(p-value) close to zero because those
  # gene sets are not particularly interesting, and there are many of
  # them, which slows down the plotting.
  i   <- which(with(dat,!(abs(p0) < 3 & abs(p1) < 3)))
  j   <- which(with(dat,abs(p0) < 3 & abs(p1) < 3))
  n   <- length(j)
  j   <- sample(j,ceiling(n/10))
  dat <- dat[c(i,j),]

  # Set any missing signed -log10(p-value) to zero.
  dat[is.na(dat$p1),"p1"] <- 0
  dat[is.na(dat$p0),"p0"] <- 0
  return(dat)
}


# Create a scatterplot to visualize the gene-set enrichment results
# for a given topic k. Input argument "gene_set_info" is a data frame
# containing information about the gene sets; "gsea_res" is an output
# from function "perform_gsea"; and "label_gene_sets" is a vector of
# gene set ids to be labeled in the plot.
create_genescore_gsea_plot <- function (gene_set_info, gsea_res, k,
                                        label_gene_sets = NULL,
                                        title = paste("topic",k)) {

  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k)

  # Remove all gene sets ids except those that were specifically
  # mentioned by the "label_gene_sets" input argument.
  if (!is.null(label_gene_sets))
    dat[!is.element(dat$id,label_gene_sets),"id"] <- NA

  # Create the scatterplot.
  return(ggplot(dat,aes_string(x = "p0",y = "p1",shape = "database",
                               fill = "collection",label = "id")) +
           geom_point(size = 2,color = "white",stroke = 0.3) +
           geom_abline(intercept = 0,slope = 1,color = "gray",
                       linetype = "dashed") +
           geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                           box.padding = 0.1,point.padding = 0.3,
                           segment.color = "black",segment.size = 0.25,
                           min.segment.length = 0.1,na.rm = TRUE) +
           scale_shape_manual(values = gsea_plot_shapes) +
           scale_fill_manual(values = gsea_plot_colors) +
           guides(shape = guide_legend(override.aes = list(fill = "black"))) +
           guides(fill = guide_legend(override.aes = list(shape = 21))) +
           labs(x = "most extreme signed -log10(p-value) among other topics",
                y = "signed -log10(p-value) in topic",
                title = title) +
           theme_cowplot(font_size = 12) +
           theme(plot.title = element_text(size = 12,face = "plain")))
}

# Create an interactive scatterplot using plotly to explore the
# gene-set enrichment results. The input arguments are similar to
# create_gsea_plot.
create_genescore_gsea_plotly <- function (gene_set_info, gsea_res, k, margin = 1,
                                          height = 675, width = 800,
                                          title = paste("topic",k)) {

  # Compile the data for plotting.
  dat <- compile_data_for_gsea_plot(gene_set_info,gsea_res,k)

  # Create the plotly plot.
  pmin <- min(c(dat$p0,dat$p1))
  pmax <- max(c(dat$p0,dat$p1))
  p <- plot_ly(data = dat,x = ~p0,y = ~p1,symbol = ~database,
               color = ~collection,
               text = ~sprintf("%s\nid: %s\np-value: %0.2e",name,id,10^(-p1)),
               type = "scatter",mode = "markers",hoverinfo = "text",
               symbols = gsea_plot_shapes,colors = gsea_plot_colors,
               marker = list(size = 8,line = list(color = "white",width = 1)),
               height = height,width = width)
  p <- add_trace(p,data = data.frame(x = c(pmin,pmax),y = c(pmin,pmax)),
                 x = ~x,y = ~y,mode = "lines",type = "scatter",
                 inherit = FALSE,showlegend = FALSE,
                 line = list(color = "lightgray",dash = "dash",size = 0.3))
  p <- layout(p,
              legend = list(font = list(size = 10)),
              xaxis = list(title="most extreme signed -log10(p-value) in other topics",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              yaxis = list(title = "signed -log10(p-value) in topic",
                           showgrid = FALSE,showline = FALSE,zeroline = FALSE),
              hoverlabel = list(bgcolor = "white",bordercolor = "black",
                                font = list(color = "black",family = "arial",
                                            size = 12)),
              font = list(family = "arial",size = 12),
              title = title)
  return(p)
}


#' Compile HOMER motif enrichment results
#'
#' @param homer_res_topics a list of Homer results,
#' each contains the Homer motif enrichment results for a topic after
#' running the `run_homer()` function in `motif_analysis.R`.
#'
#' @export
#'
compile_homer_motif_res <- function(homer_res_topics) {

  motif_names <- sort(unique(as.character(sapply(homer_res_topics, function(x) x[,1]))))

  motif_mlog10P <- matrix(NA, nrow = length(motif_names), ncol = length(homer_res_topics))
  colnames(motif_mlog10P) <- names(homer_res_topics)
  rownames(motif_mlog10P) <- motif_names

  motif_Padj <- motif_mlog10P
  colnames_homer <- c("motif_name", "consensus", "P", "logP", "Padj",  "num_target", "percent_target", "num_bg", "percent_bg")

  for (k in 1:length(homer_res_topics)){
    homer_res <- homer_res_topics[[k]]
    colnames(homer_res) <- colnames_homer
    motif_mlog10P[,k] <- -1 * homer_res[match(motif_names, homer_res$motif_name), "logP"]/log(10)
    motif_Padj[,k] <- homer_res[match(motif_names, homer_res$motif_name), "Padj"]
  }

  motif_ranks <- apply(-1*motif_mlog10P, 2, rank, ties.method = "min")

  motifs <- data.frame(x = motif_names) %>% separate(x, c("motif", "origin", "database"), "/")
  motifs$motif <- gsub("\\s*\\(.*\\)", "", motifs$motif)
  rownames(motifs) <- motif_names

  motif_logP <- -log(10)*motif_mlog10P
  motif_zscore <- pval_to_zscore(motif_logP, tails=1, log.p=TRUE)

  motif_res <- list(motifs = motifs,
                    mlog10P = motif_mlog10P,
                    Z = motif_zscore,
                    Padj = motif_Padj,
                    rank = motif_ranks)
  return(motif_res)
}


# Compile the data frame used for motif enrichment plots
compile_data_for_motif_plot <- function (motif_res, k, x, y, subsample = FALSE) {

  P <- motif_res$mlog10P
  # Compute the biggest -log10(p-value) among other topics.
  n  <- nrow(P)
  p0 <- rep(0,n)
  for (i in 1:n) {
    if (is.na(P[i,k]))
      p0[i] <- NA
    else
      p0[i] <- max(P[i,-k],na.rm = TRUE)
  }

  # Compute the lowest ranking among other topics.
  ranks <- motif_res$rank
  rank0 <- rep(0,n)
  for (i in 1:n) {
    rank0[i] <- min(ranks[i,-k],na.rm = TRUE)
  }

  motif <- gsub("\\s*\\(.*", "", motif_res$motifs$motif)

  # Compile the data for plotting.
  dat <- data.frame(motif = motif,
                    p1 = P[,k],
                    p0 = p0,
                    rank1 = ranks[,k],
                    rank0 = rank0,
                    stringsAsFactors = FALSE)

  # Subsample the motifs with -log10(p-value) close to zero because those
  # motifs are not particularly interesting.
  if (subsample) {
    i   <- which(with(dat,!(abs(p0) < 3 & abs(p1) < 3)))
    j   <- which(with(dat,abs(p0) < 3 & abs(p1) < 3))
    n   <- length(j)
    j   <- sample(j,ceiling(n/10))
    dat <- dat[c(i,j),]
  }

  return(dat)
}

# Create a scatterplot to visualize the motif enrichment results
# for a given topic k.
create_motif_enrichment_plot <- function (motif_res, k,
                                          label_motifs = NULL,
                                          title = paste("topic",k),
                                          max.overlaps = Inf, font.size = 9,
                                          subsample = FALSE) {

  # Compile the data for plotting.
  dat <- compile_data_for_motif_plot(motif_res, k, subsample)

  # Remove all gene sets ids except those that were specifically
  # mentioned by the "label_gene_sets" input argument.
  if (!is.null(label_motifs))
    dat[!is.element(dat$motif,label_motifs),"id"] <- NA

  dat$x <- dat$p0
  dat$y <- dat$p1
  x.label <- "max -log10(P-value) among other topics"
  y.label <- "-log10(P-value) in topic"

  # Create the scatterplot.
  return(ggplot(dat,aes_string(x = "x",y = "y", label = "motif")) +
           geom_point(size = 2,stroke = 0.3) +
           geom_abline(intercept = 0,slope = 1,color = "gray",
                       linetype = "dashed") +
           geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                           segment.color = "black",segment.size = 0.25,
                           max.overlaps = max.overlaps,na.rm = TRUE) +
           labs(x = x.label,
                y = y.label,
                title = title) +
           theme_cowplot(font.size))

}

# Create a scatterplot to visualize the motif enrichment results
# for a given topic k.
create_motif_enrichment_ranking_plot <- function (motif_res, k,
                                                  label_motifs = NULL,
                                                  title = paste("topic",k),
                                                  max.overlaps = Inf, font.size = 9,
                                                  subsample = FALSE) {

  # Compile the data for plotting.
  dat <- compile_data_for_motif_plot(motif_res, k, subsample)

  # Remove all gene sets ids except those that were specifically
  # mentioned by the "label_gene_sets" input argument.
  if (!is.null(label_motifs))
    dat[!is.element(dat$motif,label_motifs),"id"] <- NA

  dat$y <- dat$p1
  y.label <- "-log10(P-value)"

  # Create the scatterplot.
  return(ggplot(dat,aes_string(x = "rank1",y = "y", label = "motif")) +
           geom_point(size = 2,stroke = 0.3) +
           geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                           segment.color = "black",segment.size = 0.25,
                           max.overlaps = max.overlaps,na.rm = TRUE) +
           labs(x = "Ranking",
                y = y.label,
                title = title) +
           theme_cowplot(font.size))

}

#' Create a heatmap to visualize the motif enrichment results across topics
#'
#' @param motif_res a list of motif enrichment statistics from `compile_homer_motif_res` function.
#' @param enrichment enrichment value to plot in the heatmap, must be "z-score" or "-log10(p-value)"
#' @param cluster_motifs If true, cluster motifs by hierarchical clustering.
#' @param cluster_topics If ture, cluster topics by hierarchical clustering.
#' @param motif_filter_value Filter out motifs. Motifs inclued in the heamtmap need to have
#' maximum abs(zscore) or -log10(p-value) across topics greater than motif_filter_value.
#' @param enrichment_range Truncate values in the heamtpa to be within this range.
#' @param method_cluster method for hierarchical clustering.
#' @param horizontal If true, plot motifs on the x-axis.
#' @param color.low the low end color in the heatmap
#' @param color.high the high end color in the heatmap
#' @param font.size.motifs the font size for motif names
#' @param font.size.topics the font size for the topic names
#'
create_motif_enrichment_heatmap <- function (motif_res,
                                             enrichment = c("z-score","-log10(p-value)"),
                                             cluster_motifs = TRUE,
                                             cluster_topics = TRUE,
                                             motif_filter_value = 10,
                                             enrichment_range = c(-100,100),
                                             method_cluster = "average",
                                             horizontal = FALSE,
                                             color.low = "deepskyblue",
                                             color.high = "orangered",
                                             font.size.motifs = 6,
                                             font.size.topics = 9) {

  enrichment <- match.arg(enrichment)

  if (enrichment == "z-score"){
    enrichment.label <- "z-score"
    dat <- motif_res$Z
    rownames(dat) <- motif_res$motifs$motif
    # Filter out motifs with maximum zscore > motif_filter_value
    dat <- dat[which(apply(dat, 1, max) > motif_filter_value),]
    cat(sprintf("%d out of %d motifs included the heatmap\n", nrow(dat), nrow(motif_res$motifs)))
  }else if (enrichment == "-log10(p-value)"){
    enrichment.label <- "-log10(p-value)"
    dat <- motif_res$mlog10P
    rownames(dat) <- motif_res$motifs$motif
    # Filter out motifs with maximum -log10(p-value) > motif_filter_value
    dat <- dat[which(apply(dat, 1, max) > motif_filter_value),]
    cat(sprintf("%d out of %d motifs included the heatmap\n", nrow(dat), nrow(motif_res$motifs)))
  }

  # clustering motifs
  if (cluster_motifs) {
    motif_order <- hclust(dist(dat, method = "euclidean"), method = method_cluster)$order
  } else {
    motif_order <- order(rownames(dat))
  }

  # clustering topics
  if (cluster_topics) {
    topic_order <- hclust(dist(t(dat), method = "euclidean"), method = method_cluster)$order
  } else {
    topic_order <- 1:ncol(dat)
  }

  dat <- dat[rev(motif_order), topic_order]

  dat2 <- melt(dat)
  colnames(dat2) <- c("motif", "topic", "enrichment")

  dat2 <- transform(dat2, enrichment = pmin(max(enrichment_range), enrichment))
  dat2 <- transform(dat2, enrichment = pmax(min(enrichment_range), enrichment))

  # Heatmap
  if (horizontal) {
    p <- ggplot(dat2, aes_string(x= "motif", y = "topic", fill = "enrichment")) +
      geom_tile() +
      labs(x = "",
           y = "",
           title = "Motif enrichment",
           fill = enrichment.label) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0))
    theme(axis.text.x  = element_text(size = font.size.motifs, angle = 90),
          axis.text.y  = element_text(size = font.size.topics),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())

  }else{
    p <- ggplot(dat2, aes_string(x = "topic",y = "motif", fill = "enrichment")) +
      geom_tile() +
      labs(x = "",
           y = "",
           title = "Motif enrichment",
           fill = enrichment.label) +
      scale_x_discrete(position = "top") +
      theme(axis.text.x  = element_text(size = font.size.topics),
            axis.text.y  = element_text(size = font.size.motifs),
            legend.title=element_text(size=10),
            legend.text=element_text(size=9),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank())

  }


  if (enrichment == "z-score"){
    p <- p + scale_fill_gradient2(low = color.low,mid = "white",high = color.high, midpoint = 0)
  }else if (enrichment == "-log10(p-value)"){
    p <- p + scale_fill_gradient2(low="white", high=color.high)
  }

  return(p)

}

create_celllabel_topic_heatmap <- function (fit_multinom,
                                            grouping,
                                            group_labels,
                                            topic_order = 1:ncol(fit_multinom$L),
                                            cluster_topics = FALSE,
                                            method_cluster = "average",
                                            color.low = "white",
                                            color.mid = "orangered",
                                            color.high = "darkred",
                                            font.size.groups = 8,
                                            font.size.topics = 8) {

  L <- fit_multinom$L

  L_group_mean <- matrix(NA, nrow = length(group_labels), ncol = ncol(L))
  rownames(L_group_mean) <- group_labels
  colnames(L_group_mean) <- colnames(L)
  for (i in 1:length(group_labels)) {
    L_group <- L[grouping == group_labels[i], ]
    L_group_mean[i,] <- colMeans(L_group)
  }

  # clustering topics
  if (cluster_topics) {
    topic_order <- hclust(dist(t(L_group_mean), method = "euclidean"), method = method_cluster)$order
  } else {
    topic_order <- topic_order
  }

  dat <- L_group_mean[rev(group_labels), topic_order]

  dat2 <- melt(dat)
  colnames(dat2) <- c("group", "topic", "proportion")

  # Heatmap
  p <- ggplot(dat2, aes_string(x = "topic", y = "group", fill = "proportion")) +
    geom_tile() +
    labs(x = "",
         y = "",
         title = "",
         fill = "average proportion") +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(low = color.low, mid = color.mid, high = color.high, midpoint = 0.5) +
    theme(axis.text.x  = element_text(size = font.size.topics),
          axis.text.y  = element_text(size = font.size.groups),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())

  return(p)

}

# Create a scatterplot to visualize the motif enrichment and correlation with gene scores
# for a given topic k.
create_motif_gene_cor_scatterplot <- function (motif_res,
                                               genescore_res,
                                               motif_gene_table,
                                               k,
                                               cor.motif = c("z-score","-log10(p-value)"),
                                               cor.method = c("pearson", "spearman"),
                                               title = paste("topic",k),
                                               max.overlaps = 15,
                                               font.size = 9) {

  # Check and process input arguments.
  cor.motif <- match.arg(cor.motif)
  cor.method <- match.arg(cor.method)

  # Extract the gene scores of the selected genes
  gene_names <- genescore_res$genes$SYMBOL
  selected_genes <- motif_gene_table$gene
  idx_genes <- match(toupper(motif_gene_table$gene), toupper(gene_names))
  gene_scores_matched <- genescore_res$Z[idx_genes,]
  gene_names_matched <- gene_names[idx_genes]

  # Use the most significant motif version in this topic, if there are more than one motif version for this TF
  motif_order <- order(motif_res$mlog10P[,k], decreasing = T)
  motif_mlog10P <- motif_res$mlog10P[motif_order,]
  motif_zscore <- motif_res$Z[motif_order,]
  motif_names <- motif_res$motifs[motif_order, "motif"]

  idx_motifs <- match(toupper(motif_gene_table$motif), toupper(motif_names)) # match the first (top) match
  motif_mlog10P_matched <- motif_mlog10P[idx_motifs,]
  motif_zscore_matched <- motif_zscore[idx_motifs,]
  motif_names_matched <- gsub("/.*", "", rownames(motif_mlog10P)[idx_motifs])

  motif_gene_mapping <- data.frame(motif = motif_names_matched,
                                   gene = gene_names_matched,
                                   motif_mlog10P = motif_mlog10P_matched[,k],
                                   motif_zscore = motif_zscore_matched[,k],
                                   gene_score = gene_scores_matched[,k],
                                   cor_zscore = NA, cor_mlog10P= NA,
                                   row.names = selected_genes)

  for(i in 1:length(selected_genes)){
    motif_gene_mapping$cor_zscore[i] <- cor(motif_zscore_matched[i,], gene_scores_matched[i,], method = cor.method)
    motif_gene_mapping$cor_mlog10P[i] <- cor(motif_mlog10P_matched[i,], gene_scores_matched[i,], method = cor.method)
  }

  if(cor.motif == "z-score"){
    dat <- data.frame(x = motif_gene_mapping$cor_zscore,
                      y = motif_gene_mapping$motif_mlog10P,
                      gene = motif_gene_mapping$gene)
    x.label <- "Correlation between motif enrichment z-scores and gene scores"

  }else if (cor.motif == "-log10(p-value)"){
    dat <- data.frame(x = motif_gene_mapping$cor_mlog10P,
                      y = motif_gene_mapping$motif_mlog10P,
                      gene = motif_gene_mapping$gene)
    x.label <- "Correlation between motif enrichment -log10(p-value) and gene scores"
  }

  y.label <- "Motif enrichment -log10(p-value)"

  # Create the scatterplot.
  p <- ggplot(dat,aes_string(x = "x",y = "y", label = "gene")) +
    geom_point(size = 2,stroke = 0.3) +
    geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                    segment.color = "black",segment.size = 0.25,
                    max.overlaps = max.overlaps,na.rm = TRUE) +
    labs(x = x.label,
         y = y.label,
         title = title) +
    theme_cowplot(font.size)
  print(p)

  return(motif_gene_mapping)

}

# Create a scatterplot to visualize the motif enrichment and gene scores
# for a given topic k.
create_motif_gene_scatterplot <- function (motif_res,
                                           genescore_res,
                                           motif_gene_table,
                                           k,
                                           y = c("z-score","-log10(p-value)"),
                                           cor.method = c("pearson", "spearman"),
                                           colors,
                                           max.overlaps = 15,
                                           font.size = 9) {

  # Check and process input arguments.
  y <- match.arg(y)
  cor.method <- match.arg(cor.method)

  # Extract the gene scores of the selected genes
  gene_names <- genescore_res$genes$SYMBOL
  selected_genes <- motif_gene_table$gene
  idx_genes <- match(toupper(motif_gene_table$gene), toupper(gene_names))
  gene_scores_matched <- matrix(genescore_res$Z[idx_genes,], nrow = length(selected_genes))
  gene_names_matched <- gene_names[idx_genes]

  # Use the most significant motif version in this topic, if there are more than one motif version for this TF
  motif_order <- order(motif_res$mlog10P[,k], decreasing = T)
  motif_mlog10P <- motif_res$mlog10P[motif_order,]
  motif_zscore <- motif_res$Z[motif_order,]
  motif_names <- motif_res$motifs[motif_order, "motif"]
  idx_motifs <- match(toupper(motif_gene_table$motif), toupper(motif_names)) # match the first (top) match

  motif_mlog10P_matched <- matrix(motif_mlog10P[idx_motifs,], nrow = length(selected_genes))
  motif_zscore_matched <- matrix(motif_zscore[idx_motifs,], nrow = length(selected_genes))
  motif_names_matched <- gsub("/.*", "", rownames(motif_mlog10P)[idx_motifs])

  plots <- vector("list", length(selected_genes))
  names(plots) <- selected_genes

  for( i in 1:length(selected_genes)){
    gene.label <- gene_names_matched[i]

    if (y == "-log10(p-value)"){
      dat <- data.frame(x = gene_scores_matched[i,],
                        y = motif_mlog10P_matched[i,],
                        topics = factor(colnames(genescore_res$Z), levels = colnames(genescore_res$Z)))
      y.label <- "motif enrichment -log10(p-value)"
      title <- sprintf("%s (r = %.2f)", gene.label, cor(dat$x, dat$y, method = cor.method))
    }else if (y == "z-score"){
      dat <- data.frame(x = gene_scores_matched[i,],
                        y = motif_zscore_matched[i,],
                        topics = factor(colnames(genescore_res$Z), levels = colnames(genescore_res$Z)))
      y.label <- "motif enrichment z-score"
      title <- sprintf("%s (r = %.2f)", gene.label, cor(dat$x, dat$y, method = cor.method))
    }

    x.label <- "gene score"

    plots[[i]] <- ggplot(dat,aes_string(x = "x",y = "y", label = "topics", color = "topics")) +
      geom_point(size = 2,stroke = 0.3) +
      scale_color_manual(values = colors) +
      geom_text_repel(color = "black",size = 2.25,fontface = "italic",
                      segment.color = "black",segment.size = 0.25,
                      max.overlaps = max.overlaps,na.rm = TRUE) +
      labs(x = x.label,
           y = y.label,
           title = title) +
      theme_cowplot(font.size)
  }

  return(plots)
}

# convert p-value or logP to z-score
pval_to_zscore <- function(p, direction=NULL, tails=2, log.p=FALSE, p.limit=.Machine$double.xmin) {

  if (log.p == TRUE){
    logP <- p
  }else{
    if ( !is.null( p.limit ) ){
      p[which(p < p.limit )] <- p.limit ## set lower limit to avoid Inf/-Inf zscores
    }
    logP <- log(p)
  }

  if (tails == 2) {
    z <- qnorm(logP - log(2), lower.tail = FALSE, log.p = TRUE)
  } else if (tails == 1){

    if ( !is.null( p.limit ) ){
      logP[which(logP == 0)] <- -p.limit ## avoid -Inf zscores when logP = 0 (pvalue = 1)
    }

    z <- qnorm(logP, lower.tail = FALSE, log.p = TRUE)
  } else {
    stop( "Parameter 'tails' must be set to either 1 or 2.")
  }

  if ( !is.null( direction) ) {
    z <-  z * sign( direction )
  }

  return(z)
}

# Make motif logo plots using `Logolas` package
plot_motif_logo <- function(homer_res_topics, motif_name, k, motif.dir,
                            type = c("Logo", "EDLogo", "both"),
                            colors = c("green", "blue", "orange", "red"),
                            xaxis_fontsize = 6, y_fontsize = 6, main_fontsize=8, ...) {

  type <- match.arg(type)

  homer_res <- homer_res_topics[[k]]
  colnames(homer_res) <- c("motif_name", "consensus", "P", "logP", "Padj",  "num_target", "percent_target", "num_bg", "percent_bg")

  idx_motif <- which(homer_res$motif_name == motif_name)
  motif_file <- paste0(motif.dir, "/known", idx_motif, ".motif")
  if (!file.exists(motif_file)) {
    cat(sprintf("The PWM of the motif (%s) was not in HOMER output of enriched motifs.\n", motif_name))
    return(NULL)
  }
  motif_pwm <- read.table(motif_file, header = F, comment.char = ">")
  colnames(motif_pwm) <- c("A", "C", "G", "T")

  motif_pwm <- t(motif_pwm)

  motif <- gsub("/.*", "", motif_name)
  pval  <- homer_res[idx_motif, "P"]
  mlog10P <- -1*homer_res[idx_motif, "logP"]/log(10)
  # title <- sprintf("%s -log10(p-value) = %.1f", motif, mlog10P)

  if (type == "Logo") {
    logomaker(motif_pwm, type = "Logo", color_type = "per_row", colors = colors,
              logo_control = list(pop_name = paste(motif, "logo"),
                                  xaxis_fontsize = xaxis_fontsize, y_fontsize = y_fontsize, main_fontsize = main_fontsize))
  }else if (type == "EDLogo") {
    logomaker(motif_pwm, type = "EDLogo", color_type = "per_row", colors = colors,
              logo_control = list(pop_name = paste(motif, "ED logo"),
                                  xaxis_fontsize = xaxis_fontsize, y_fontsize = y_fontsize, main_fontsize = main_fontsize))
  }else if (type == "both") {
    Logolas::get_viewport_logo(1, 2)
    seekViewport(paste0("plotlogo", 1))
    logomaker(motif_pwm, type = "Logo", color_type = "per_row", colors = colors,
              logo_control = list(newpage = FALSE, pop_name = paste(motif, "logo"),
                                  xaxis_fontsize = xaxis_fontsize, y_fontsize = y_fontsize, main_fontsize = main_fontsize))

    seekViewport(paste0("plotlogo", 2))
    logomaker(motif_pwm, type = "EDLogo", color_type = "per_row", colors = colors,
              logo_control = list(newpage = FALSE, pop_name = paste(motif, "ED logo"),
                                  xaxis_fontsize = xaxis_fontsize, y_fontsize = y_fontsize, main_fontsize = main_fontsize))
  }

}


#' Compute the gene scores per topic/cluster using the "TSS model"
#' Weight chromatin accessibility around gene promoters by distance to TSS.
#' Adopted from Lareau et al. Nature Biotech, 2019.
#'
#' @param Z A matrix containing the scores (e.g. Z-scores) for the accessibility regions in each topic or cluster.
#' Rows are ATAC-seq regions, columns are topics or clusters.
#' @param regions A data frame containing the coordinates (chr, start, end) of the accessibility regions.
#' @param genes A data frame containing gene coordinates (chr, start, end, strand, gene_name, etc. ) .
#' @param normalize logical indicating if the weights of the regions match to each gene should be normalized.
#' @param method.normalization Normalization method. "l2": normalize by the l2 norm of the weights (default).
#' "sum": normalize by the sum of weights.
#' @param c scaling constant (default = 5000)
#' @param gene.window size of the gene window around TSS (default = 100000, i.e. 100kb)
#'
#' @return Returns a matrix of gene scores. Rows are genes, columns are topics or clusters.
#' @export
compute_gene_scores_tss_model <- function(Z,
                                          regions,
                                          genes,
                                          normalize = TRUE,
                                          method.normalization = "l2",
                                          c = 5000,
                                          window.size = 100000
) {

  # Get ATAC-seq regions
  regions <- as.data.frame(regions)
  if(nrow(Z) != nrow(regions)){
    stop("The number of regions should match with the number of rows in Z!")
  }

  colnames(regions)[1:3] <- c("chr", "start", "end")
  regions <- regions %>% mutate_at(c("start", "end"), as.numeric)
  cat(sprintf("load %d regions. \n", nrow(regions)))

  # Represent the region/peaks with the region/peak centers
  regions$center <- (regions$start + regions$end)/2
  regions <- makeGRangesFromDataFrame(regions, start.field = "center", end.field = "center")

  # Get gene windows around TSS
  genes <- as.data.frame(genes)
  colnames(genes)[1:5] <- c("chr", "start", "end", "strand", "gene_name")
  genes <- genes %>% mutate_at(c("start", "end"), as.numeric)
  cat(sprintf("load %d genes \n", nrow(genes)))

  genes$TSS <- NA
  genes$TSS[which(genes$strand == "+")] <- genes$start[which(genes$strand == "+")]
  genes$TSS[which(genes$strand == "-")] <- genes$end[which(genes$strand == "-")]

  gene.windows <- genes
  gene.windows$start <- ifelse(genes$TSS - window.size > 0, genes$TSS - window.size, 0)
  gene.windows$end <- genes$TSS + window.size

  gene.windows <- makeGRangesFromDataFrame(gene.windows, keep.extra.columns = TRUE)

  cat("Compute gene scores. \n")

  # Overlap between ATAC-seq regions and genes
  overlaps <- as.data.frame(findOverlaps(gene.windows, regions))
  colnames(overlaps) <- c("gene", "region")

  # Compute weights using the simple distance decay model from TSS
  distTSS <- mcols(gene.windows)$TSS[overlaps$gene] - start(regions)[overlaps$region]
  overlaps$weight <- exp(-abs(distTSS/c))

  # Weight matrix for genes x regions
  W <- Matrix::sparseMatrix(i = overlaps$gene,
                            j = overlaps$region,
                            x = overlaps$weight,
                            dims = c(length(gene.windows), length(regions)))
  rownames(W) <- mcols(gene.windows)$gene_name

  # Compute gene scores

  # filter out genes that do not match to any regions
  W <- W[which(Matrix::rowSums(W) != 0), ]

  # set gene scores as the weighted sum
  Z.genescore <- W %*% Z

  # normalized scores
  if (normalize == TRUE) {
    if (method.normalization == "sum") {
      # normalize by the sum of weights
      Z.genescore <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% Z.genescore
    } else {
      # normalize by the l2 norm of weights, as in Stouffer's z-score method
      Z.genescore <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(W^2))) %*% Z.genescore
    }
  }

  return(Z.genescore)

}

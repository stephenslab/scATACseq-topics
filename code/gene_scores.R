
#' Compute the gene scores per topic/cluster using the "TSS model"
#' Weight chromatin accessibility around gene promoters by distance to TSS.
#' Adopted from Lareau et al. Nature Biotech, 2019.
#'
#' @param region.scores A matrix containing the scores (e.g. Z-scores) for the accessibility regions in each topic or cluster.
#' Rows are ATAC-seq regions, columns are topics or clusters.
#' @param regions A data frame containing the coordinates (chr, start, end) of the accessibility regions.
#' @param genes A data frame or `GRanges` object containing gene coordinates (chr, start, end, strand, gene_name, etc. ) .
#' @param c scaling constant (default = 5000)
#' @param gene.window size of the gene window around TSS (default = 100000, i.e. 100kb)
#'
#' @return Returns a matrix of gene scores. Rows are genes, columns are topics or clusters.
#' @export
compute_gene_scores_tss_model <- function(region.scores,
                                          regions,
                                          genes,
                                          c = 5000,
                                          window.size = 100000
) {

  # Get ATAC-seq regions/peaks
  regions <- as.data.frame(regions)
  cat(sprintf("load %d regions. \n", nrow(regions)))

  if(nrow(region.scores) != nrow(regions)){
    stop("The numbers of regions do not match!")
  }

  colnames(regions)[1:3] <- c("chr", "start", "end")
  regions <- regions %>% mutate_at(c("start", "end"), as.numeric)

  # Use ATAC-seq region/peak centers to represent the region/peaks
  regions$center <- (regions$start + regions$end)/2
  regions <- makeGRangesFromDataFrame(regions, start.field = "center", end.field = "center")

  # Get gene windows around TSS
  genes <- as.data.frame(genes)
  colnames(genes)[1:3] <- c("chr", "start", "end")
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
  match_region2gene <- as.data.frame(findOverlaps(regions, gene.windows))
  colnames(match_region2gene) <- c("region", "gene")

  # Compute weights using the simple distance decay model from TSS
  distTSS <- mcols(gene.windows)$TSS[match_region2gene$gene] - start(regions)[match_region2gene$region]

  match_region2gene$weight <- exp(-abs(distTSS/c))

  # Weight matrix for regions x genes
  W <- Matrix::sparseMatrix(i = match_region2gene$region,
                            j = match_region2gene$gene,
                            x = match_region2gene$weight,
                            dims = c(length(regions), length(gene.windows)))
  colnames(W) <- genes$gene_name

  # Compute gene scores with the weights
  gene.scores <- t(W) %*% region.scores

  return(gene.scores)

}

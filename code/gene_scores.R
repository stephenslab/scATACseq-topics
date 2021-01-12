
#' Compute the gene scores per topic/cluster using the "TSS model"
#' Weight chromatin accessibility around gene promoters by distance to TSS.
#' The "TSS model" uses bi-directional exponential decays from the gene TSS.
#' Adopted from Lareau et al. Nature Biotech, 2019.
#'
#' @param Z A matrix containing the scores (e.g. Z-scores) for the accessibility regions in each topic or cluster.
#' Rows are ATAC-seq regions, columns are topics or clusters.
#' @param ATAC.regions A data frame containing the coordinates (chr, start, end) of the accessibility regions.
#' @param genes A data frame containing gene coordinates (chr, start, end, strand, GeneID, etc. ).
#' @param use.ATAC.centers logical indicating whether to represent ATAC-seq positions by the centers of ATAC-seq regions
#' @param normalize logical indicating if the weights of the regions match to each gene should be normalized.
#' @param method.normalization Normalization method (`l2` or `sum`). `l2: normalize by the l2 norm of the weights (default).
#' `sum`: normalize by the sum of weights.
#' @param c scaling constant (default = 5000)
#' @param window.upstream An integer specifying the size of the window upstream of TSS (default = 100000, i.e. 100kb)
#' @param window.downstream An integer specifying the size of the window downstream of TSS (default = 100000, i.e. 100kb)
#'
#' @return Returns a matrix of gene scores. Rows are genes, columns are topics or clusters.
#' @export
compute_gene_scores_tss_model <- function(Z,
                                          ATAC.regions,
                                          genes,
                                          use.ATAC.centers = TRUE,
                                          normalize = TRUE,
                                          method.normalization = "l2",
                                          c = 5000,
                                          window.upstream = 100000,
                                          window.downstream = 100000
) {

  if(is.matrix(Z)){
    Z <- as.matrix(Z)
  }

  # Get ATAC-seq regions
  ATAC.regions <- as.data.frame(ATAC.regions)
  if(nrow(Z) != nrow(ATAC.regions)){
    stop("The number of ATAC regions should match with the number of rows in Z!")
  }

  colnames(ATAC.regions)[1:3] <- c("chr", "start", "end")
  ATAC.regions <- ATAC.regions %>% mutate_at(c("start", "end"), as.numeric)

  if (use.ATAC.centers) {
    # Represent the region/peaks with the region/peak centers
    ATAC.regions$center <- (ATAC.regions$start + ATAC.regions$end)/2
    ATAC.regions <- makeGRangesFromDataFrame(ATAC.regions, start.field = "center", end.field = "center")
  }else {
    ATAC.regions <- makeGRangesFromDataFrame(ATAC.regions, start.field = "start", end.field = "end")
  }

  # Get gene windows around TSS
  genes <- as.data.frame(genes)
  colnames(genes)[1:5] <- c("chr", "start", "end", "strand", "GeneID")
  genes <- genes %>% mutate_at(c("start", "end"), as.numeric)

  genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

  # Get gene TSS
  genes.TSS <- resize(genes, 1, "start")

  # Get gene windows by extending upstream (and downstream) around the TSS
  gene.windows <- extend_ranges(genes.TSS, window.upstream, window.downstream)

  # Find all ATAC-seq regions within the gene window
  cat(sprintf("Match %d genes with %d ATAC regions. \n", length(gene.windows), length(ATAC.regions)))
  overlaps <- as.data.frame(findOverlaps(gene.windows, ATAC.regions))
  colnames(overlaps) <- c("idx_gene", "idx_ATAC")

  cat("Compute gene scores. \n")

  # Compute distance from ATAC-seq regions to gene TSS
  dist <- distance_ranges(ATAC.regions[overlaps$idx_ATAC], genes.TSS[overlaps$idx_gene])

  # Compute weight matrix for genes x ATAC regions
  weight <- exp(-abs(dist/c))

  W <- Matrix::sparseMatrix(i = overlaps$idx_gene,
                            j = overlaps$idx_ATAC,
                            x = weight,
                            dims = c(length(gene.windows), length(ATAC.regions)))
  rownames(W) <- mcols(gene.windows)$GeneID

  # Filter out genes that do not match to any ATAC regions
  W <- W[which(Matrix::rowSums(W) != 0), ]

  # Compute the weighted sum of region scores
  Z.genescore <- W %*% Z

  # Normalize gene scores
  if (normalize == TRUE) {
    if (method.normalization == "l2") {
      # normalize by the l2 norm of weights, as in Stouffer's z-score method
      Z.genescore <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(W^2))) %*% Z.genescore
    } else {
      # normalize by the sum of weights
      Z.genescore <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% Z.genescore
    }
  }

  return(Z.genescore)

}




#' Compute the gene scores per topic/cluster using the gene score model (model 42)
#' based on the archR paper with some modifications.
#' modified based on https://github.com/GreenleafLab/ArchR/blob/master/R/MatrixGeneScores.R
#' The gene score model uses bi-directional exponential decays from the gene TSS (extended upstream by 5 kb by default)
#' and the gene transcription termination site (TTS).
#' Note: the current version of the function does not account for neighboring gene boundaries,
#' and does not perform scaling by gene size.
#'
#' @param Z A matrix containing the scores (e.g. Z-scores) for the accessibility regions in each topic or cluster.
#' Rows are ATAC-seq regions, columns are topics or clusters.
#' @param ATAC.regions A data frame containing the coordinates (chr, start, end) of the accessibility regions.
#' @param genes A data frame containing gene coordinates (chr, start, end, strand, GeneID, etc. ) .
#' @param use.ATAC.centers logical indicating whether to represent ATAC-seq positions by the centers of ATAC-seq regions
#' @param weight.model A string for the weighting model for weighting ATAC-seq regions for gene score calculation.
#' This string should be a function of `dist`, where `dist` is the distance from the ATAC-seq regions to the gene.
#' Default gene model: "exp(-abs(dist)/5000) + exp(-1)". dist is the distance to gene body.
#' @param distTo A string, genebody (default) or TSS. `genebody` will compute distances from the ATAC-seq regions to gene body.
#' `TSS` will compute distances from the ATAC-seq regions to TSS.
#' @param normalize logical indicating if the weights of the regions match to each gene should be normalized.
#' @param method.normalization Normalization method (`l2` or `sum`). `l2: normalize by the l2 norm of the weights (default).
#' `sum`: normalize by the sum of weights.
#' @param window.upstream An integer specifying the size of the window upstream of TSS (default = 100000, i.e. 100kb)
#' @param window.downstream An integer specifying the size of the window downstream of TSS (default = 100000, i.e. 100kb)
#' @param gene.upstream An integer describing the number of bp upstream the gene to extend the gene body (default = 5000).
#' @param gene.downstream An integer describing the number of bp downstream the gene to extend the gene body (default = 0).
#'
#' @return Returns a matrix of gene scores. Rows are genes, columns are topics or clusters.
#' @export
compute_gene_scores_genebody_model <- function(Z,
                                               ATAC.regions,
                                               genes,
                                               use.ATAC.centers = TRUE,
                                               weight.model = "exp(-abs(dist)/5000) + exp(-1)",
                                               distTo = "genebody",
                                               normalize = TRUE,
                                               method.normalization = "l2",
                                               window.upstream = 100000,
                                               window.downstream = 100000,
                                               gene.upstream = 5000,
                                               gene.downstream = 0
) {

  if(is.matrix(Z)){
    Z <- as.matrix(Z)
  }

  # Get ATAC-seq regions
  ATAC.regions <- as.data.frame(ATAC.regions)
  if(nrow(Z) != nrow(ATAC.regions)){
    stop("The number of ATAC regions should match with the number of rows in Z!")
  }

  colnames(ATAC.regions)[1:3] <- c("chr", "start", "end")
  ATAC.regions <- ATAC.regions %>% mutate_at(c("start", "end"), as.numeric)

  if (use.ATAC.centers) {
    # Represent the region/peaks with the region/peak centers
    ATAC.regions$center <- (ATAC.regions$start + ATAC.regions$end)/2
    ATAC.regions <- makeGRangesFromDataFrame(ATAC.regions, start.field = "center", end.field = "center")
  }else {
    ATAC.regions <- makeGRangesFromDataFrame(ATAC.regions, start.field = "start", end.field = "end")
  }

  # Get gene windows around TSS
  genes <- as.data.frame(genes)
  colnames(genes)[1:5] <- c("chr", "start", "end", "strand", "GeneID")
  genes <- genes %>% mutate_at(c("start", "end"), as.numeric)

  genes <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

  # Get gene TSS
  genes.TSS <- resize(genes, 1, "start")

  # Get expanded gene body by extending upstream (and downstream) around the genes
  genes.genebody <- extend_ranges(genes, gene.upstream, gene.downstream)

  # Get gene windows by extending upstream (and downstream) around the TSS
  gene.windows <- extend_ranges(genes.TSS, window.upstream, window.downstream)

  # Find all ATAC-seq regions within the gene window
  cat(sprintf("Match %d genes with %d ATAC regions. \n", length(gene.windows), length(ATAC.regions)))

  overlaps <- as.data.frame(findOverlaps(gene.windows, ATAC.regions))
  colnames(overlaps) <- c("idx_gene", "idx_ATAC")

  # Compute distance from ATAC-seq regions to genes
  cat("Compute gene scores. \n")
  if (distTo == "TSS") {
    dist <- distance_ranges(ATAC.regions[overlaps$idx_ATAC], genes.TSS[overlaps$idx_gene])
  }else{
    dist <- distance_ranges(ATAC.regions[overlaps$idx_ATAC], genes.genebody[overlaps$idx_gene])
  }

  # Compute weight matrix for genes x ATAC regions
  weight <- eval(parse(text=weight.model))

  W <- Matrix::sparseMatrix(i = overlaps$idx_gene,
                            j = overlaps$idx_ATAC,
                            x = weight,
                            dims = c(length(gene.windows), length(ATAC.regions)))
  rownames(W) <- mcols(gene.windows)$GeneID

  # Filter out genes that do not match to any ATAC regions
  W <- W[which(Matrix::rowSums(W) != 0), ]

  # Compute the weighted sum of region scores
  Z.genescore <- W %*% Z

  # Normalize gene scores
  if (normalize == TRUE) {
    if (method.normalization == "l2") {
      # normalize by the l2 norm of weights, as in Stouffer's z-score method
      Z.genescore <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(W^2))) %*% Z.genescore
    } else {
      # normalize by the sum of weights
      Z.genescore <- Matrix::Diagonal(x = 1 / Matrix::rowSums(W)) %*% Z.genescore
    }
  }

  return(Z.genescore)

}

#' Extend upstream and downstream from a GenomicRanges object
#'
#' @param x A GenomicRanges object.
#' @param upstream The number of nucleotides toward the 5' end of each region.
#' @param downstream The number of nucleotides toward the 3' end of each region.
#' @export
extend_ranges <- function(x, upstream=2000, downstream=2000, ignore.strand = FALSE) {

  if (ignore.strand == FALSE) {
    is_minus_strand <- which(strand(x) == "-")
    is_plus_strand <- which(strand(x) != "-")

    start(x)[is_plus_strand] <- start(x)[is_plus_strand] - upstream
    end(x)[is_plus_strand] <- end(x)[is_plus_strand] + downstream

    end(x)[is_minus_strand] <- end(x)[is_minus_strand] + upstream
    start(x)[is_minus_strand] <- start(x)[is_minus_strand] - downstream
  }else{
    start(x) <- start(x) - upstream
    end(x) <- end(x) + downstream
  }

  start(x)[which(start(x) < 1)] <- 1

  return(x)
}

#' Compute the genomic distance between ATAC-seq regions and genes
#'
#' @param regions A GenomicRanges object of the ATAC-seq regions or other genomic loci.
#' @param genes A GenomicRanges object of the genes
#' @export
distance_ranges <- function(regions, genes) {
  absDist <- distance(regions, genes)
  signDist <- sign(start(regions) - start(resize(genes,1,"start")))

  is_minus_strand <- which(strand(genes) == "-")
  signDist[is_minus_strand] <- signDist[is_minus_strand] * -1
  dist <- absDist * signDist

  return(dist)
}

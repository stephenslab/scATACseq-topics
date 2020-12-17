# Some unctions for Gene set enrichment analysis
# borrowed from https://github.com/stephenslab/single-cell-topics/blob/master/code/gsea.R

# This function aligns the gene-set data (gene_sets) with the
# gene-wise statistics (gene_scores) by their gene IDs to
# prepare these data for a gene-set enrichment analysis.
# Modified based on the `align_gene_data` function in `single-cell-topics` repo.

# It is assumed that the row names of gene_sets and the row names of
# gene_scores matrices give the unique gene symbols or ids (e.g., Ensembl ids).
align_gene_data <- function (gene_sets, gene_scores) {
  x   <- rownames(gene_sets)
  y   <- rownames(gene_scores)
  ids <- intersect(x,y)
  i   <- match(ids,x)
  j   <- match(ids,y)
  gene_sets   <- gene_sets[i,]
  gene_scores <- gene_scores[j,]
  return(list(gene_sets = gene_sets, gene_scores = gene_scores))
}

# Recover a list of gene sets from n x m adjacency matrix A, in which
# n is the number of genes and m is the number of gene sets, and
# A[i,j] is greater than zero if and only if gene i is included in
# gene set j. This function can be used to prepare the "pathways"
# input to fgsea from a collection of gene sets encoded as a matrix.
# Note that the rows and columns of the matrix, A, should be named.
matrix2pathways <- function (A) {
  n          <- ncol(A)
  out        <- vector("list",n)
  names(out) <- colnames(A)
  genes      <- rownames(A)
  for (i in 1:n)
    out[[i]] <- genes[A[,i] > 0]
  return(out)
}

# Perform gene-set enrichment analysis using the fast method described
# in Korotkevich et al (2016). This method uses permutation testing
# to account for correlations between genes. The p-values are based on
# Kolmogorov-Smirnov test statistics which are "normalized" to account
# for differences in gene-set sizes; the normalized test statistics
# are provided in the "NES" (normalized enrichment score) column of
# the output. See also Subramanian et al (2005) for background on the
# method.
#
# Input argument "gene_sets" should be an n x m matrix, in which n is
# the number of genes, and m is the number of gene sets, and entry
# (i,j) is greater than zero if and only if gene i is included in gene
# set j. Input argument z should be a numeric vector of length n
# containing gene-wise statistics, such as z-scores. The elements of
# z should be named.
#
# Input argument "eps" controls the accuracy of the small p-values.
# Here I set it to be much lower than the default setting since some
# of the gene-set enrichment p-values can be quite small.
perform_gsea <- function (gene_sets, z, eps = 1e-32, nproc = 1, ...) {

  # Convert the gene-sets adjacency matrix into the fgsea gene-sets
  # format.
  pathways <- matrix2pathways(gene_sets)

  # Perform the gene-set enrichment analysis using fgsea.
  out <- suppressWarnings(fgsea(pathways,z,eps = eps,nproc = nproc,...))
  class(out) <- "data.frame"

  # Post-process the fgsea output.
  rownames(out) <- out$pathway
  out <- out[c("pval","log2err","ES","NES")]
  out <- out[colnames(gene_sets),]
  out[is.na(out$log2err),] <- NA
  return(out)
}

# Perform fgsea gene-set enrichment analysis once per column of the
# gene_scores matrix.
# Modified based on the `perform_gsea_all_topics` function in `single-cell-topics` repo to
# use Z instead of diff_count_res as input argument.
#' @param gene_sets n x m matrix, in which n is
# the number of genes, and m is the number of gene sets, and entry
# (i,j) is greater than zero if and only if gene i is included in gene
# set j.
#' @param Z n x k matrix of gene statistics (e.g. z-scores).
#' n is the number of genes and k is the number of topics.
#' @export
perform_gsea_all_topics <- function (gene_sets, Z,
                                     eps = 1e-32, nproc = 1, ...) {

  # Get the number of gene sets (n) and the number of topics (k).
  n <- ncol(gene_sets)
  k <- ncol(Z)

  # Initialize the outputs.
  out <- list(pval    = matrix(0,n,k),
              log2err = matrix(0,n,k),
              ES      = matrix(0,n,k),
              NES     = matrix(0,n,k))
  rownames(out$pval)    <- colnames(gene_sets)
  rownames(out$log2err) <- colnames(gene_sets)
  rownames(out$ES)      <- colnames(gene_sets)
  rownames(out$NES)     <- colnames(gene_sets)
  colnames(out$pval)    <- colnames(Z)
  colnames(out$log2err) <- colnames(Z)
  colnames(out$ES)      <- colnames(Z)
  colnames(out$NES)     <- colnames(Z)

  # Run the gene-set enrichment analysis for each topic.
  cat("k = ")
  for (i in 1:k) {
    if (i == 1)
      cat("1")
    else
      cat(",",i)
    z               <- Z[,i]
    names(z)        <- rownames(Z)
    ans             <- perform_gsea(gene_sets,z,eps,nproc,...)
    out$pval[,i]    <- ans$pval
    out$log2err[,i] <- ans$log2err
    out$ES[,i]      <- ans$ES
    out$NES[,i]     <- ans$NES
  }
  cat("\n")

  # Output the results of the gene-set enrichment analysis.
  return(out)
}

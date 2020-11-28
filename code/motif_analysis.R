## Run HOMER for motif enrichment

#' Finding Enriched Motifs in Genomic Regions using HOMER
#'
#' @param path.to.homer Directory path to 'findMotifsGenome.pl' excutable file.
#' @param region.bed BED file of target regions
#' @param genome Genome version (default='hg19').
#' @param result.dir
#' @param region.size The size of the region used for motif finding (default=200)
#' @param motif.length Specifies the length of motifs to be found (default=8,10,12).
#' HOMER will find motifs of each size separately and then combine the results at the end.
#' The length of time it takes to find motifs increases greatly with increasing size.
#' @param optimize.count Mismatches allowed in global optimization phase (default=2)
#' @param n.motifs Number of motifs to find (default=25)
#' @param use.hypergeometric use hypergeometric distribution to score motif enrichment
#' By default, findMotifsGenome.pl uses the binomial distribution to score motifs, which is faster.
#' This works well when the number of background sequences greatly out number the target sequences.
#' If the number of background sequences is smaller than target sequences, it is recommended to use the hypergeometric distribution.
#' @param background.bed BED file of background regions
#' @param n.cores Number of cores to use.
#'
#' @return
#' @export
#'
#' @examples
homer_motif_enrichment <- function(path.to.homer = 'findMotifsGenome.pl',
                                   region.bed,
                                   genome = 'hg19',
                                   result.dir,
                                   region.size = 200,
                                   motif.length = c(8,10,12),
                                   optimize.count = 2,
                                   n.motifs = 25,
                                   use.hypergeometric = FALSE,
                                   background.bed = NULL,
                                   n.cores = 5){

  path.to.homer <- normalizePath(path.to.homer)

  if (!file_test('-x', path.to.homer)) {
    stop(path.to.homer, " does not exist or is not executable; check your path.to.homer parameter")
  }

  ## Make HOMER results output dir
  dir.create(result.dir, showWarnings = FALSE, recursive = TRUE)

  cmd <- paste(
    path.to.homer,
    region.bed, genome, result.dir,
    '-len', paste0(motif.length, collapse = ','),
    '-size', region.size,
    '-mis', optimize.count,
    '-S', n.motifs,
    '-p', n.cores
    )

  if (use.hypergeometric == TRUE) {
    cmd <- paste(cmd, '-h')
  }

  if (!is.null(background.bed)) {
    cmd <- paste(cmd, '-bg', background.bed)
  }

  cat("Run HOMER command: \n", cmd)

  system(cmd)
  x <- read.csv(paste0(result.dir, "/knownResults.txt"), sep="\t", header=TRUE)
  system("rm -f *.tmp")
  return(x)
}

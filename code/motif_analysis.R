
#' Run HOMER for motif enrichment
#' See http://homer.ucsd.edu/homer/ngs/peakMotifs.html for more details.
#' @param path.to.homer Directory path to 'findMotifsGenome.pl' excutable file.
#' @param regions BED file of the peak regions.
#' @param genome Genome assembly (default='hg19').
#' @param out.dir Output directory
#' @param region.size Size of the regions used for motif finding (default=200)
#' @param motif.length Length of motifs to be found (default=8,10,12).
#' HOMER will find motifs of each size separately and then combine the results at the end.
#' The length of time it takes to find motifs increases greatly with increasing size.
#' @param optimize.count Number of mismatches allowed in global optimization phase (default=2)
#' @param n.motifs Number of motifs to find (default=25)
#' @param use.hypergeometric Whether to use hypergeometric distribution to score motif enrichment
#' By default, findMotifsGenome.pl uses the binomial distribution to score motifs, which is faster.
#' This works well when the number of background sequences greatly out number the target sequences.
#' If the number of background sequences is smaller than target sequences, it is recommended to use the hypergeometric distribution.
#' @param background BED file of the background regions
#' @param n.cores Number of cores to use.
#' @param run Whether to run this command inside R (default=TRUE)
#' @export
run_homer <- function(path.to.homer = 'findMotifsGenome.pl',
                      regions = "peaks.bed",
                      genome = 'hg19',
                      out.dir = "out",
                      region.size = 200,
                      motif.length = c(8,10,12),
                      optimize.count = 2,
                      n.motifs = 25,
                      use.hypergeometric = FALSE,
                      background = NULL,
                      n.cores = 1,
                      run = TRUE
) {

  path.to.homer <- normalizePath(path.to.homer)

  if(!file_test('-x', path.to.homer)){
    stop(path.to.homer, " does not exist or is not executable; check your path.to.homer parameter")
  }

  ## Make HOMER results output dir
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

  cmd <- paste(
    path.to.homer,
    regions, genome, out.dir,
    '-len', paste0(motif.length, collapse = ','),
    '-size', region.size,
    '-mis', optimize.count,
    '-S', n.motifs,
    '-p', n.cores
  )

  if(use.hypergeometric){
    cmd <- paste(cmd, '-h')
  }

  if(!is.null(background)){
    cmd <- paste(cmd, '-bg', background)
  }

  cat("\nHOMER command: \n")
  cat(cmd, "\n")

  if(run){
    system(cmd)
    x <- read.csv(file.pah(out.dir, "knownResults.txt"), sep="\t", header=TRUE)
    system("rm -f *.tmp")
    return(x)
  }

}

#' Run GREAT for genomic regions ontologies enrichment
#' For more details, see https://jokergoo.github.io/rGREAT and https://bioconductor.org/packages/release/bioc/html/rGREAT.html
#' @param gr A GRanges object or a data frame which contains at least three columns (chr, start and end). Regions for test.
#' @param bg A GRanges object or a data frame. Background regions if needed. Note gr should be exactly subset of bg for all columns in gr.
#' @param genome Genome assembly. "hg38", "hg19", "mm10", "mm9" are supported in GREAT version 4.x.x,
#' "hg19", "mm10", "mm9", "danRer7" are supported in GREAT version 3.x.x.
#' @param includeCuratedRegDoms Whether to include curated regulatory domains.
#' @param rule How to associate genomic regions to genes: "basalPlusExt", "twoClosest" or "oneClosest".
#' @param request_interval Time interval for two requests.
#' @param max_tries Maximum times trying to connect to GREAT web server.
#' @param version version of GREAT (Default = 4.0)
#' @param base_url the url of cgi-bin path, only used when explicitly specified.
#' @export
run_GREAT <- function(gr,
                      bg = NULL,
                      genome,
                      includeCuratedRegDoms = TRUE,
                      rule = "basalPlusExt",
                      request_interval = 300,
                      max_tries = 10,
                      version = 4.0,
                      base_url = "http://great.stanford.edu/public/cgi-bin"
) {
  library(rGREAT)
  job = submitGreatJob(gr,
                       bg                    = bg,
                       species               = genome,
                       includeCuratedRegDoms = includeCuratedRegDoms,
                       rule                  = rule,
                       adv_upstream          = 5.0,
                       adv_downstream        = 1.0,
                       adv_span              = 1000.0,
                       adv_twoDistance       = 1000.0,
                       adv_oneDistance       = 1000.0,
                       request_interval = request_interval,
                       max_tries = max_tries,
                       version = version,
                       base_url = base_url
  );

  job
  tb = getEnrichmentTables(job)
  names(tb)

  return(list(job = job, tb = tb))
}

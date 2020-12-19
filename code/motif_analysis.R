
#' Run HOMER for motif enrichment
#' See http://homer.ucsd.edu/homer/ngs/peakMotifs.html for more details.
#' @param regions.file BED file of the peak regions.
#' @param genome Genome assembly (default='hg19').
#' @param homer.path Directory path to 'findMotifsGenome.pl' excutable file.
#' @param out.dir Output directory.
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
run_homer <- function(regions.file = "peaks.bed",
                      genome = 'hg19',
                      homer.path = 'findMotifsGenome.pl',
                      out.dir = "out",
                      region.size = 200,
                      motif.length = c(8,10,12),
                      optimize.count = 2,
                      n.motifs = 25,
                      use.hypergeometric = FALSE,
                      background = NULL,
                      n.cores = 1,
                      dryrun = FALSE
) {

  cmd <- paste(
    homer.path,
    regions.file, genome, out.dir,
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

  if(dryrun){
    cat(sprintf("HOMER command:\n%s\n", cmd))
  }else{
    if(!dir.exists(out.dir))
      dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
    homer.path <- normalizePath(homer.path)
    if(!file_test('-x', homer.path)){
      stop(homer.path, " does not exist or is not executable!")
    }
    system(cmd)
    res <- read.csv(file.path(out.dir, "knownResults.txt"), sep="\t", header=TRUE)
    return(res)
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


#' Select regions based on differential accessbility result and save selected regions as BED files.
#' @param diff_count_res differential accessbility result from `diff_count_analysis`
#' @param method method to select regions, `FDR` selects based on FDR, `logFC` selects based on logFC,
#' `FDR&logFC` selects based on FDR and logFC, quantile selectes regions above abs(z-score) quantile.
#'
#' @param our.dir Output directory.
#' @param thresh.FDR Threshold of FDR, default = 0.01.
#' @param thresh.logFC Threshold of abs(logFC), default = 1.
#' @param thresh.quantile Threshold of abs(z-score) quantile, default = 0.99.
#' @param save.bed If TRUE, save selected regions as BED files.
select_regions <- function(diff_count_res,
                           method="FDR",
                           out.dir = "out",
                           thresh.FDR = 0.01,
                           thresh.logFC = 1,
                           thresh.quantile = 0.99,
                           save.bed = TRUE) {

  selected_regions <- vector("list", length = ncol(diff_count_res$Z))
  names(selected_regions) <- colnames(diff_count_res$Z)
  selected_regions$filenames <- c()
  for(k in colnames(diff_count_res$Z)){
    z <- diff_count_res$Z[,k]
    beta <- diff_count_res$beta[,k]
    nlog10P <- diff_count_res$pval[,k]
    p <- 10^(-nlog10P)

    if(toupper(method) == "FDR"){
      idx_regions_sig <- which(p.adjust(p, method = "fdr") < thresh.FDR)
      cat(sprintf("%s: %d regions with FDR < %.3f. \n", k, length(idx_regions_sig), thresh.FDR))
    }else if(tolower(method) == "bonferroni"){
      idx_regions_sig <- which(p.adjust(p, method = "bonferroni") < thresh.FDR)
      cat(sprintf("%s: %d regions with Bonferroni p.adj < %.3f. \n", k, length(idx_regions_sig), thresh.FDR))
    }else if(method == "qvalue"){
      library(qvalue)
      idx_regions_sig <- which(qvalue(p)$qvalues < thresh.FDR)
      cat(sprintf("%s: %d regions with qvalue < %.3f. \n", k, length(idx_regions_sig), thresh.FDR))
    }else if(method == "logFC"){
      idx_regions_sig <- which(abs(beta) > thresh.logFC)
      cat(sprintf("%s: %d regions with abs(logFC) > %.1f. \n", k, length(idx_regions_sig), thresh.logFC))
    }else if(method == "quantile"){
      idx_regions_sig <- which(abs(z) > quantile(abs(z), thresh.quantile))
      cat(sprintf("%s: %d regions above quantile of %.3f. \n", k, length(idx_regions_sig), thresh.quantile))
    }else if(method == "FDR_logFC"){
      idx_regions_sig <- which( (p.adjust(p, method = "BH") < thresh.FDR) & (abs(beta) > thresh.logFC) )
      cat(sprintf("%s: %d regions with FDR < %.3f AND abs(logFC) > %.3f.\n", k, length(idx_regions_sig), thresh.FDR, thresh.logFC))
    }else{
      stop("Method not recognized!")
    }

    if(length(idx_regions_sig) > 0){
      region_names <- names(z)[idx_regions_sig]
      regions <- data.frame(x = region_names)
      regions <- regions %>% separate(x, c("chr", "start", "end"), "_")  %>% mutate_at(c("start", "end"), as.numeric)
      regions <- data.frame(regions, name = region_names)
      selected_regions[[k]] <- regions
      if(save.bed){
        # save selected marker regions as BED files
        if(!dir.exists(out.dir))
          dir.create(out.dir, showWarnings = FALSE, recursive = T)
        bedfile <- paste0(out.dir, "/selected_regions_", k, ".bed")
        write.table(regions, file = bedfile, quote=F, sep="\t", row.names=F, col.names=F)
        selected_regions$filenames <- c(selected_regions$filenames, bedfile)
      }
    }

  }

  return(selected_regions)
}

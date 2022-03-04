
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
#' @return A data frame of motif enrichment results loaed from HOMER's output 'knownResults.txt'
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
    res <- read.csv(file.path(out.dir, "knownResults.txt"), sep="\t", header=TRUE, check.names=F, stringsAsFactors=F)
    return(res)
  }

}

#' Select regions based on differential accessibility result and save selected regions as BED files.
#' @param DA_res Differential accessibility result from `de_analysis`.
#' @param method Method to select regions.
#' `pval` selects regions in which p-value < `thresh.pval`.
#' `lfsr` selects regions in which lfsr > `thresh.lfsr`.
#' `logFC` selects regions in which beta > `thresh.logFC`.
#' `topPercent` selects top regions in higest `topPercent`,
#' `topN` selects the top `n.regions` regions.
#' @param out.dir Output directory.
#' @param thresh.pval p-value threshold.
#' @param thresh.lfsr lfsr threshold.
#' @param thresh.logFC logFC threshold.
#' @param top.percent Select top percent regions with largest logFC.
#' @param top.n Select top N regions with largest logFC.
#' @param save.bed If TRUE, save selected regions as BED files for downstream motif enrichment analysis.
#' @import dplyr tidyr
#' @return A list of selected regions for each topic
#' @export
select_DA_regions <- function(DA_res,
                              method = c("pval", "lfsr", "logFC", "topPercent", "topN"),
                              out.dir = "out",
                              thresh.pval = 0.1,
                              thresh.lfsr = 0.1,
                              thresh.logFC = 2,
                              top.percent = 0.01,
                              top.n = 2000,
                              save.bed = TRUE) {

  method <- match.arg(method)

  selected_regions <- vector("list", length = ncol(DA_res$z))
  names(selected_regions) <- colnames(DA_res$z)
  selected_regions$filenames <- c()
  for(k in colnames(DA_res$z)){
    z <- na.omit(DA_res$z[,k])
    logFC <- na.omit(DA_res$postmean[,k])
    lpval <- na.omit(DA_res$lpval[,k])
    lfsr <- na.omit(DA_res$lfsr[,k])

    if(method == "pval"){
      sig_regions <- which(lpval > -log10(thresh.pval))
      cat(sprintf("topic %s: %d regions selected (pval < %.2f) \n", k, length(sig_regions), thresh.pval))
    }else if(method == "lfsr"){
      sig_regions <- which(lfsr < thresh.lfsr)
      cat(sprintf("topic %s: %d regions selected (lfsr < %.2f) \n", k, length(sig_regions), thresh.lfsr))
    }else if(method == "logFC"){
      sig_regions <- which(logFC > thresh.logFC)
      cat(sprintf("topic %s: %d regions selected by logFC > %.2f \n", k, length(sig_regions), thresh.logFC))
    }else if(method == "topPercent"){
      sig_regions <- which(logFC > quantile(logFC, 1-top.percent))
      cat(sprintf("topic %s: %d regions selected (top %.2f %%) \n", k, length(sig_regions), top.percent*100))
    }else if(method == "topN"){
      sig_regions <- head(order(logFC, decreasing = T), top.n)
      cat(sprintf("topic %s: %d regions selected by top %d regions. \n", k, length(sig_regions), top.n))
    }else{
      stop("Method not recognized! Please check 'method'.")
    }

    if(length(sig_regions) > 0){
      region_names <- names(z)[sig_regions]
      regions <- data.frame(x = region_names)
      regions <- regions %>%
        tidyr::separate(x, c("chr", "start", "end"), "_") %>%
        dplyr::mutate_at(c("start", "end"), as.numeric)
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

#' Test motif enrichment
#'
#' @param targetValue Number of successes in target
#' @param numTargets Total number in target
#' @param bgValue Number of successes in background
#' @param numBackground Total number in background
#' @param method The statistical test to perform,
#' must be one of "binomial", "hypergeometric" or "normal"
#'
#' @export
test_motif_enrichment <- function(targetValue, numTargets,
                                  bgValue, numBackground,
                                  method = c("binomial", "hypergeometric", "normal")){
  method <- match.arg(method)

  cat(sprintf("%.2f%% in target",targetValue/numTargets*100))
  cat(sprintf("%.2f%% in background",bgValue/numBackground*100))

  bgProb <- bgValue/numBackground

  if(method == "binomial"){
    logP <- pbinom(targetValue - 1, numTargets, bgProb, lower.tail = FALSE, log.p = TRUE)
  }else if(method == "hypergeometric"){
    logP <- phyper(targetValue - 1, # Number of Successes the -1 is due to cdf integration
                   bgValue, # Number of all successes in background
                   numBackground - bgValue, # Number of non successes in background
                   numTargets, # Number that were drawn
                   lower.tail = FALSE, log.p = TRUE)
  }else if(method == "normal"){
    mu <- numTargets*bgProb
    sigma <- sqrt(numTargets*bgProb*(1-bgProb))
    z <- (targetValue-mu)/sigma
    logP <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
  }

  return(logP)
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
#' @import rGREAT
#' @export
run_great <- function(gr,
                      bg = NULL,
                      genome,
                      includeCuratedRegDoms = TRUE,
                      rule = "basalPlusExt",
                      request_interval = 300,
                      max_tries = 10,
                      version = 4.0,
                      base_url = "http://great.stanford.edu/public/cgi-bin"
) {
  job <- submitGreatJob(gr,
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
  tb <- getEnrichmentTables(job)
  names(tb)

  return(list(job = job, tb = tb))
}


#' Create a gene annotation
#' Modified based on the `createGeneAnnotation` function in `ArchR`.
#'
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9"). If `genome` is not supplied,
#' `TxDb` and `OrgDb` are required. If genome is supplied, `TxDb` and `OrgDb` will be ignored.
#' @param TxDb A `TxDb` object (transcript database) from Bioconductor which contains information for gene/transcript coordinates.
#' For example, from `txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene`.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids.
#' For example, from `orgdb <- org.Hs.eg.db`.
#' @param annoStyle annotation style to map between gene names and various gene identifiers e.g. "ENTREZID", "ENSEMBL".
#' @export
get.genes.from.TxDb <- function(
  TxDb = NULL,
  OrgDb = NULL,
  annoStyle = NULL
){

  ###########################
  message("Getting Genes from TxDb...")
  genes <- GenomicFeatures::genes(TxDb)

  ## Get gene symbol from OrgDb
  if(is.null(annoStyle)){
    isEntrez <- mcols(genes)$symbol <- tryCatch({
      suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
      TRUE
    }, error = function(x){
      FALSE
    })

    isEnsembl <- mcols(genes)$symbol <- tryCatch({
      suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
      TRUE
    }, error = function(x){
      FALSE
    })

    if(isEntrez){
      annoStyle <- "ENTREZID"
    }else if(isEnsembl){
      annoStyle <- "ENSEMBL"
    }else{
      stop("Could not identify keytype for annotation format!")
    }

  }
  annoStyle <- toupper(annoStyle)

  message("Determined Annotation Style = ", annoStyle)

  mcols(genes)$symbol <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = annoStyle, multiVals = "first"))
  mcols(genes)$symbol[is.na(mcols(genes)$symbol)] <- paste0("NA_", mcols(genes)$gene_id)[is.na(mcols(genes)$symbol)]
  names(genes) <- NULL
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

  return(genes)

}

#' Get TxDb transcript metadata for the specified genome from Bioconductor.
#' Modified based on the `.getTxDb` function in `ArchR`.
#'
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9").
#' @param install Install the bioconductor TxDb package for the specified genome.
.getTxDb <- function(genome = NULL, install = TRUE){

  if(toupper(genome) == "HG19"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg19.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg19.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if(toupper(genome) == "HG38"){
    if(suppressWarnings(!require(TxDb.Hsapiens.UCSC.hg38.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", update=FALSE)
      }else{
        stop("TxDb.Hsapiens.UCSC.hg38.knownGene is not installed!")
      }
    }
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }else if(toupper(genome) == "MM9"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm9.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm9.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  }else if(toupper(genome) == "MM10"){
    if(suppressWarnings(!require(TxDb.Mmusculus.UCSC.mm10.knownGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", update=FALSE)
      }else{
        stop("TxDb.Mmusculus.UCSC.mm10.knownGene is not installed!")
      }
    }
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if(toupper(genome) == "SACCER3"){
    if(suppressWarnings(!require(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", update=FALSE)
      }else{
        stop("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene is not installed!")
      }
    }
    library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
    txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
  }else if(toupper(genome) == "RHEMAC8"){
    if(suppressWarnings(!require(TxDb.Mmulatta.UCSC.rheMac8.refGene))){
      if(install){
        message("Package does not exist, now trying bioconductor..")
        BiocManager::install("TxDb.Mmulatta.UCSC.rheMac8.refGene", update=FALSE)
      }else{
        stop("TxDb.Mmulatta.UCSC.rheMac8.refGene is not installed!")
      }
    }
    library(TxDb.Mmulatta.UCSC.rheMac8.refGene)
    txdb <- TxDb.Mmulatta.UCSC.rheMac8.refGene
  }else{
    stop("Genome not recognized!")
  }

  return(txdb)

}

#' Get OrgDb to map Entrez Gene identifiers and GenBank accession numbers.
#' Modified based on the `.getOrgDb` function in `ArchR`.
#'
#' @param genome A string that specifies the genome (ie "hg38", "hg19", "mm10", "mm9").
.getOrgDb <- function(genome = NULL){

  if(toupper(genome) == "HG19" | toupper(genome) == "HG38"){
    if(suppressWarnings(!require(org.Hs.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Hs.eg.db", update=FALSE)
    }
    library(org.Hs.eg.db)
    annodb <- org.Hs.eg.db
  }else if(toupper(genome) == "MM9" | toupper(genome) == "MM10"){
    if(suppressWarnings(!require(org.Mm.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Mm.eg.db", update=FALSE)
    }
    library(org.Mm.eg.db)
    annodb <- org.Mm.eg.db
  }else{
    stop("Genome not recognized!")
  }

  return(annodb)

}


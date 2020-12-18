#' Create a gene annotation
#' Modified based on the `createGeneAnnotation` function in `ArchR`.

#' @param TxDb A `TxDb` object (transcript database) from Bioconductor which contains information for gene/transcript coordinates.
#' For example, from `txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene`.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids.
#' For example, from `orgdb <- org.Hs.eg.db`.
#' @param keytype_geneID Keytype of the gene ID from OrgDb (e.g. "ENTREZID", "ENSEMBL").
#' @param columns_extract Other gene info columns to extract (e.g. "SYMBOL", "ENSEMBL").
#' @export
get.gene.annotations <- function(
  TxDb = NULL,
  OrgDb = NULL,
  keytype_geneID = "ENTREZID",
  columns_extract = c("SYMBOL", "ENSEMBL")
){

  ###########################
  cat("Get genes from TxDb...")
  genes        <- GenomicFeatures::genes(TxDb)
  geneID_map   <- map.geneIDs(OrgDb, genes$gene_id, keytype_geneID, columns_extract)
  mcols(genes) <- cbind(mcols(genes), geneID_map[,columns_extract])
  names(genes) <- NULL
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

  return(genes)

}


#' Extract mapped gene IDs from OrgDb
#' Modified based on the `createGeneAnnotation` function in `ArchR`.
#' @param OrgDb An `OrgDb` object (organism database) from Bioconductor which contains information for gene/transcript symbols from ids.
#' For example, from `orgdb <- org.Hs.eg.db`.
#' @param geneIDs Input gene IDs.
#' @param keytype_geneID Keytype of the input gene IDs (e.g. "ENTREZID", "ENSEMBL").
#' @param columns_extract Other gene ID columns to extract (e.g. "SYMBOL", "ENSEMBL").
#' @export
map.geneIDs <- function(
  OrgDb = NULL,
  geneIDs = NULL,
  keytype_geneID = NULL,
  columns_extract = c("SYMBOL", "ENSEMBL")
){

  # Determine the keytype of gene IDs if not provided
  if(is.null(keytype_geneID)){
    isEntrez <- tryCatch({
      suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
      TRUE
    }, error = function(x){
      FALSE
    })

    isEnsembl <- tryCatch({
      suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = geneIDs, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
      TRUE
    }, error = function(x){
      FALSE
    })

    if(isEntrez){
      keytype_geneID <- "ENTREZID"
    }else if(isEnsembl){
      keytype_geneID <- "ENSEMBL"
    }else{
      stop("Could not identify keytype for annotation format!")
    }
  }
  keytype_geneID <- toupper(keytype_geneID)

  cat(sprintf("Input keytype of the gene IDs: %s \n", keytype_geneID))

  # Extract the other types of gene ID columns mapped to the input gene IDs
  columns_extract <- setdiff(toupper(columns_extract), keytype_geneID)
  geneID_map <- data.frame(matrix(nrow = length(geneIDs), ncol = length(columns_extract)))
  colnames(geneID_map) <- columns_extract

  for(column in columns_extract){
    cat(sprintf("Extract: %s \n", column))
    geneID_map[,column] <- suppressMessages(AnnotationDbi::mapIds(OrgDb, keys = geneIDs, column = column, keytype = keytype_geneID, multiVals = "first"))
    geneID_map[is.na(geneID_map[,column]), column] <- paste0("NA_", geneIDs[is.na(geneID_map[,column])])
  }

  return(geneID_map)

}

# load TxDb for human or mouse
# modified based on https://github.com/GreenleafLab/ArchR/blob/master/R/AnnotationGenome.R
getTxDb <- function(genome = NULL, install = TRUE){
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
  }else{
    stop("Genome not recognized!")
  }

  return(txdb)

}

# load OrgDb for human or mouse
# modified based on https://github.com/GreenleafLab/ArchR/blob/master/R/AnnotationGenome.R
getOrgDb <- function(genome = NULL){
  if(toupper(genome) == "HG19" | toupper(genome) == "HG38"){
    if(suppressWarnings(!require(org.Hs.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Hs.eg.db", update=FALSE)
    }
    library(org.Hs.eg.db)
    OrgDb <- org.Hs.eg.db
  }else if(toupper(genome) == "MM9" | toupper(genome) == "MM10"){
    if(suppressWarnings(!require(org.Mm.eg.db))){
      message("Package does not exist, now trying bioconductor..")
      BiocManager::install("org.Mm.eg.db", update=FALSE)
    }
    library(org.Mm.eg.db)
    OrgDb <- org.Mm.eg.db
  }else{
    stop("Genome not recognized!")
  }
  return(OrgDb)

}


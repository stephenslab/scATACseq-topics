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

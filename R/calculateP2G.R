#' Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param ArchR_path String specifying the path to an ArchR project if ArchR's implementation of addPeak2GeneLinks is desired
#' @param useDim String specifying the dimensional reduction representation in the ArchR project to use or the name of the reduced dimension matrix supplied by the user
#' @param useMatrix String specifying which the name of the gene expression matrix in the ArchR project to use.
#' It is often called the "GeneExpressionMatrix" for multiome and "GeneIntegrationMatrix" for unpaired data in ArchR project.
#' @param cor_cutoff A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
#'  Default correlation cutoff is 0.5.
#' @param ... other parameters to pass to addPeak2GeneLinks from ArchR package
#'
#' @return A DataFrame of Peak to Gene correlation
#' @import SummarizedExperiment SingleCellExperiment GenomicRanges
#' @export
#'
#' @author Xiaosai Yao, Shang-yang Chen

calculateP2G <- function(ArchR_path = NULL,
                         useDim = "IterativeLSI",
                         useMatrix = "GeneIntegrationMatrix",
                         cor_cutoff = 0.5,
                         ...) {


  if (!is.null(ArchR_path)) {


    writeLines("Using ArchR to compute peak to gene links...")
    suppressMessages(obj <- ArchR::loadArchRProject(ArchR_path))

    obj <- ArchR::addPeak2GeneLinks(
      ArchRProj = obj,
      reducedDims = useDim,
      useMatrix = useMatrix,
      logFile = "x",
      ...
    )

    p2g <- ArchR::getPeak2GeneLinks(
      ArchRProj = obj,
      corCutOff = cor_cutoff,
      resolution = 1000,
      returnLoops = FALSE
    )

    # Get metadata from p2g object and turn into df with peak indexes
    peak_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[1]]) # shows  chromosome, start, and end coordinates for each peak
    peak_metadata$idxATAC <- seq_along(rownames(peak_metadata))

    gene_metadata <- as.data.frame(S4Vectors::metadata(p2g)[[2]]) # shows gene name and RNA index of genomic ranges
    gene_metadata$idxRNA <- seq_along(rownames(gene_metadata))

    # Add gene names and peak positions to dataframe
    p2g_merged <- merge(p2g, gene_metadata, by = "idxRNA") # merge by gene ID
    p2g_merged <- merge(p2g_merged, peak_metadata, by = "idxATAC") # merge by peak ID

    # Calculate distance
    p2g_merged$distance <- abs((p2g_merged$start.y + p2g_merged$end.y)/2 - (p2g_merged$start.x + p2g_merged$end.x)/2)

    # Extract relevant columns
    p2g_merged <- p2g_merged[, c("idxATAC", "seqnames.y", "start.y","end.y", "idxRNA", "name", "Correlation", "distance")]
    colnames(p2g_merged) <- c("idxATAC", "chr", "start","end", "idxRNA", "target", "Correlation", "distance")




  } else {
    stop(
      "Input obj must be an 'ArchR' path"
    )

  }
  p2g_merged <- p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA),,drop=FALSE]
  return(p2g_merged)

}




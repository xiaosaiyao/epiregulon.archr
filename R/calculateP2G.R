#' Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
#'
#' @param peakMatrix A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq.
#' `rowRanges` should contain genomic positions of the peaks in the form of `GRanges`.
#' @param expMatrix A SingleCellExperiment object containing gene expression counts from scRNA-seq. `rowRanges` should contain genomic positions of
#' the genes in the form of `GRanges`. `rowData` should contain a column of gene symbols with column name matching the `gene_symbol` argument.
#' @param reducedDim A matrix of dimension reduced values
#' @param ArchR_path String specifying the path to an ArchR project if ArchR's implementation of addPeak2GeneLinks is desired
#' @param cor_cutoff A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
#'  Default correlation cutoff is 0.5.
#' @param useDim String specifying the dimensional reduction representation in the ArchR project to use or the name of the reduced dimension matrix supplied by the user
#' @param useMatrix String specifying which the name of the gene expression matrix in the ArchR project to use.
#' It is often called the "GeneExpressionMatrix" for multiome and "GeneIntegrationMatrix" for unpaired data in ArchR project.
#' @param ... other parameters to pass to [ArchR::addPeak2GeneLinks] package or to [epiregulon::calculateP2G].
#'
#' @return A DataFrame of Peak to Gene correlation
#' @details Cluster information is sometimes helpful to avoid the [Simpsons's paradox](https://en.wikipedia.org/wiki/Simpson%27s_paradox) in which baseline differences
#' between cell lines or cell types can create artificial or even inverse correlations between peak accessibility and gene expression. If Cluster information is provided,
#' correlation is performed within cell aggregates of each cluster.
#' @import SummarizedExperiment SingleCellExperiment GenomicRanges
#' @export
#'
#' @examples
#' # create a mock singleCellExperiment object for gene expression matrix
#' set.seed(1000)
#' gene_sce <- scuttle::mockSCE()
#' gene_sce <- scuttle::logNormCounts(gene_sce)
#' gene_gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3","chr4"), nrow(gene_sce)/4),
#'                    ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
#'                    width = 100))
#' rownames(gene_sce) <- rownames(gene_sce)
#' gene_gr$name <- rownames(gene_sce)
#' rowRanges(gene_sce) <- gene_gr
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
#'                       nrow = length(peak_gr), ncol=ncol(gene_sce))
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak",1:10)

#' # create a mock reducedDim matrix
#' reducedDim_mat <- matrix(runif(ncol(gene_sce)*50, min = 0, max = 1), nrow = ncol(gene_sce), 50)
#' p2g <- calculateP2G(peakMatrix = peak_sce, expMatrix = gene_sce, reducedDim = reducedDim_mat,
#'                     cellNum = 20, clusters = gene_sce$Treatment)
#' @author Xiaosai Yao, Shang-yang Chen

calculateP2G <- function(peakMatrix = NULL,
                         expMatrix = NULL,
                         reducedDim = NULL,
                         ArchR_path = NULL,
                         useDim = "IterativeLSI",
                         useMatrix = "GeneIntegrationMatrix",
                         cor_cutoff = 0.5,
                         ...) {


  if (!is.null(ArchR_path)) {
    ArchR::addArchRLogging(useLogs = FALSE)

    writeLines("Using ArchR to compute peak to gene links...")
    suppressMessages(obj <- ArchR::loadArchRProject(ArchR_path))
    additional_arguments <- list(...)[c(names(list(...)) %in% names(formals(ArchR::getPeak2GeneLinks)))]
    obj <- do.call(ArchR::addPeak2GeneLinks, c(list(ArchRProj = obj, reducedDims = useDim,
                                                       useMatrix = useMatrix, logFile = "x"),
                                               additional_arguments))

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
    p2g_merged <- p2g_merged[order(p2g_merged$idxATAC, p2g_merged$idxRNA),,drop=FALSE]
    return(p2g_merged)

  } else if (!is.null(peakMatrix) &
             !is.null(expMatrix) & !is.null(reducedDim)) {
    additional_arguments <- list(...)[c(names(list(...)) %in% names(formals(epiregulon::calculateP2G)))]
    return(do.call(epiregulon::calculateP2G, c(list(peakMatrix = peakMatrix, expMatrix = expMatrix,
                                                       reducedDim = reducedDim, useDim = useDim,
                                                       cor_cutoff = cor_cutoff, BPPARAM = BiocParallel::SerialParam()),
                                               additional_arguments)))

  } else {
    stop(
      "Input obj must be either an 'ArchR' path or all 3 matrices: gene expression, chromatin accessibility and dimensionality reduction"
    )

  }

}

#' Retrieve TF binding sites or motif positions
#'
#' Combined transcription factor ChIP-seq data from ChIP-Atlas and ENCODE or
#' from CistromeDB and ENCODE.
#' @param mode a string indicating whether to download a GRangelist of TF binding sites ("occupancy") or motif matches ("motif").
#' TF binding information is retrieved from [scMultiome](https://github.com/xiaosaiyao/scMultiome/blob/devel/R/tfBinding.R) whereas
#' motif information is annotated by cisbp from [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs)
#' @param archr_path Character string indicating the path of the ArchR project to retrieve motif information if
#' motif enrichment was already performed. If no motif enrichment has been performed, first annotate the ArchR using
#' `addMotifAnnotations`. If no ArchR project is provided, the user can also provide peaks in the form of GRanges and
#' this function will annotate the peaks with Cisbp
#' @param motif_name Character string	indicating name of the peakAnnotation object (i.e. Motifs) to retrieve from the designated ArchRProject.
#' @param peaks A GRanges object indicating the peaks to perform motif annotation on if ArchR project is not provided.
#' The peak indices should match the `re` column in the regulon

#' @inherit scMultiome::tfBinding params return references
#' @examples
#' # retrieve TF binding info
#' \dontrun{
#' getTFMotifInfo("mm10", "atlas")
#' }
#'
#' # retrieve motif info
#' peaks <- GRanges(seqnames = c("chr12","chr19","chr19","chr11","chr6"),
#' ranges = IRanges(start = c(124914563,50850844, 50850844, 101034172, 151616327),
#' end = c(124914662,50850929, 50850929, 101034277, 151616394)))
#' grl <- getTFMotifInfo(genome = "hg38", mode = "motif", peaks=peaks)
#'
#' @export
#'
getTFMotifInfo <- function (genome = c("hg38", "hg19", "mm10"),
                            source = c("atlas", "cistrome", "encode.sample", "atlas.sample","atlas.tissue"),
                            metadata = FALSE,
                            mode = c("occupancy","motif"),
                            archr_path = NULL,
                            motif_name = "Motif",
                            peaks = NULL) {
  genome <- match.arg(genome)
  source <- match.arg(source)
  mode <- match.arg(mode)
  if (all(c(mode == "motif", is.null(peaks), is.null(archr_path)))) stop("'peaks' should be provided if 'motif' mode is chosen")
  if (mode == "occupancy") {
    grl <- scMultiome::tfBinding(genome, source, metadata)
  } else if (mode == "motif") {
    if (!is.null(archr_path)) {
      message("retrieving previously annotated motif from ArchR project")
      ArchProj <-
        ArchR::loadArchRProject(path = archr_path, showLogo = FALSE)
      matches <- ArchR::getMatches(ArchProj, name = motif_name)
      motifs <- assay(matches, "matches")

      grl <- list()
      for (tf in colnames(motifs)){
        grl[[tf]] <- SummarizedExperiment::rowRanges(matches[motifs[,tf],])
      }
      grl <- GenomicRanges::GRangesList(grl)
      names(grl) <- unlist(lapply(strsplit(names(grl), split="_|\\."), "[", 1))

    } else {
      species <- switch(genome,
                        hg38 = "human",
                        hg19 = "human",
                        mm10 = "mouse")
      BS.genome <- switch(genome,
                          hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                          hg19 = "BSgenome.Hsapiens.UCSC.hg19",
                          mm10 = "BSgenome.Mmusculus.UCSC.mm10")

      message("keeping only standard chromosomes..")
      peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")

      message("annotating peaks with motifs")
      grl <- epiregulon:::annotateMotif(species, peaks, BS.genome, out = "positions")
      names(grl) <- lapply(strsplit(names(grl), split="_|\\."), "[", 3)
    }

  }

  grl
}



#' Add TF binding motif occupancy information to the peak2gene object
#'
#' @param p2g A Peak2Gene dataframe created by ArchR or getP2Glinks() function
#' @param grl GRangeList object containing reference TF binding information. We recommend retrieving `grl` from `getTFMotifInfo` which
#' contains TF occupancy data derived from public and ENCODE ChIP-seq peaks. Alternatively, if the users would like to provide a GRangeList
#' of motif annotations. This can be derived using `motifmatchr::matchMotifs`. See details
#' @param peakMatrix A matrix of scATAC-seq peak regions with peak ids as rows
#' @param archR_project_path Path to an ArchR project that have performed LSI dimensionality reduction and scRNA-seq integration
#'
#' @return A data frame containing overlapping ids of scATAC-seq peak regions and reference TF binding regions
#' @details This function annotates each regulatory element with possible transcription factors. We can either provide a GRangeList of known ChIP-seq
#' binding sites (TF occupancy) or positions of TF motifs (TF motifs). While public ChIP-seq data may not fully align with the ground truth TF occupancy in users' data
#' (due to technical challenges of ChIP-seq or cell type nature of TF occupancy), it does offer a few important advantages over TF motif information:
#' \enumerate{
#' \item{TF occupancy allows co-activators to be included. Co-activators are chromatin modifiers that do not directly bind to DNA but nonetheless play an important role
#' in gene regulation}
#' \item{TF occupancy can distinguish between members of the same class that may share similar motifs but that may have drastically different binding sites}
#' }
#' If multiple ChIP-seq are available for the same TF, we merge the ChIP-seq data to represent an universal set of possible binding sites. The predicted TF
#' occupancy is further refined by \code{\link{pruneRegulon}}.
#'
#' If the users prefer to use TF motifs instead of TF occupancy, the users can create a GRangeList of motif annotation using `motifmatchr::matchMotifs`.
#' Here, we demonstrate how to annotate peaks with cisbp motif database
#' ```
#' library(motifmatchr)
#' library(chromVARmotifs)
#' data("human_pwms_v1")
#' peaks <- GenomicRanges::GRanges(seqnames = c("chr1","chr2","chr2"),
#'                                ranges = IRanges::IRanges(start = c(76585873,42772928, 100183786),
#'                                                          width = 500))
#' grl <- matchMotifs(human_pwms_v1, peaks, genome = "hg38", out = "positions")
#' # retain only TF symbols. TF symbols need to be consistent with gene names in regulon
#' names(grl) <- sapply(strsplit(names(grl), "_"), "[",3)
#' ```
#'
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' # create a mock peak-to-gene matrix
#' p2g <- data.frame(idxATAC = c(rep(1,5), rep(2,5)), Chrom = "chr1", idxRNA = 1:10,
#' Gene = paste0("Gene_",1:10),Correlation = runif(10, 0,1))
#'
#' # create mock a GRanges list of TF binding sites
#' grl <- GRangesList("TF1" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(50,1050), width = 100)),
#' "TF2" = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(1050), width = 100))
#' )
#'
#' # create a mock singleCellExperiment object for peak matrix
#' peak_gr <- GRanges(seqnames = "chr1",
#'              ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000),
#'              width = 100))
#' peak_counts <- matrix(sample(x = 0:4, size = 100*length(peak_gr), replace = TRUE),
#' nrow = length(peak_gr), ncol = 100)
#' peak_sce <- SingleCellExperiment(list(counts = peak_counts))
#' rowRanges(peak_sce) <- peak_gr
#' rownames(peak_sce) <- paste0("peak",1:10)
#'
#' # create overlaps between p2g matrix, TF binding sites and peak matrix
#' overlap <- addTFMotifInfo(p2g, grl, peakMatrix = peak_sce)
#' utils::head(overlap)
#' @author Xiaosai Yao, Shang-yang Chen

addTFMotifInfo <- function(p2g,
                           grl,
                           peakMatrix = NULL,
                           archR_project_path = NULL){

  if (!is.null(archR_project_path)) {
    proj <- ArchR::loadArchRProject(path = archR_project_path, showLogo = FALSE)
    peakSet <- ArchR::getPeakSet(ArchRProj = proj)
  } else {
    peakSet <- rowRanges(peakMatrix)
  }

  message("Computing overlap...")
  overlap <- GenomicRanges::findOverlaps(peakSet, grl)
  overlap <- data.frame(overlap)
  colnames(overlap) <- c("idxATAC", "idxTF")
  overlap <- overlap[which(overlap$idxATAC %in% p2g$idxATAC), ]
  overlap$tf <- names(grl)[overlap$idxTF]
  message("Success!")

  return(overlap)

}





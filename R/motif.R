#' Add Motif Scores
#'
#' @param regulon A DataFrame consisting of tf (regulator) and target in the column names.
#' @param archr_path Character string indicating the path of the ArchR project to retrieve motif information if
#' motif enrichment was already performed. If no motif enrichment has been performed, first annotate the ArchR using
#' `addMotifAnnotations`. If no ArchR project is provided, the user can also provide peaks in the form of GRanges and
#' this function will annotate the peaks with Cisbp
#' @param ArchProj An ArchR project as an alternative to providing an ArchR path
#' @param field_name Character string	indicating the column name of the regulon to add the motif information to
#' @param motif_name Character string	indicating name of the peakAnnotation object (i.e. Motifs) to retrieve from the designated ArchRProject.
#' @param peaks A GRanges object indicating the peaks to perform motif annotation on if ArchR project is not provided.
#' The peak indices should match the `re` column in the regulon
#' @param pwms A PWMatrixList for annotation of motifs using 'motifmatchr::matchMotifs'
#' @param species Character string indicating species. Currently supported species is human or mouse
#' @param genome Character string indicating the genomic build
#' @param ... Additional arguments to pass into motifmatchr::matchMotifs
#'
#' @return A DataFrame with motif matches added with 1s indicating the presence of motifs and
#' 0s indicating the absence of motifs
#' @export
#'
#' @examples
#' regulon <- S4Vectors::DataFrame(tf = c("AR","AR","AR","ESR1","ESR1","NKX2-1"),
#' idxATAC = 1:6)

#' peaks <- GRanges(seqnames = c("chr12","chr19","chr19","chr11","chr6","chr1"),
#' ranges = IRanges(start = c(124914563,50850845, 50850844, 101034172, 151616327, 1000),
#' end = c(124914662,50850929, 50850929, 101034277, 151616394,2000)))
#' regulon <- addMotifScore(regulon, peaks=peaks)



addMotifScore <- function(regulon,
                          archr_path=NULL,
                          ArchProj=NULL,
                          field_name="motif",
                          motif_name="Motif",
                          peaks=NULL,
                          pwms=NULL,
                          species=c("human","mouse"),
                          genome=c("hg38", "hg19", "mm10"),
                          ...) {

  species <- match.arg(species)
  genome <- match.arg(genome)

  if (!is.null(archr_path)){
    ArchProj <-
      ArchR::loadArchRProject(path = archr_path, showLogo = FALSE)
  }

  if (!is.null(ArchProj) & is.null(peaks)) {
    message("retrieving motif information from ArchR project")
    matches <- ArchR::getMatches(ArchProj, name = motif_name)
    motifs <- assay(matches, "matches")

    # Convert motifs to gene names
    motif_names <- unlist(lapply(strsplit(colnames(motifs), split="_"), "[", 1))

    # match motif_names and official gene symbols

    colnames(motifs) <- epiregulon:::matchNames(motif_names, regulon)
    peaks.idx <- seq_len(nrow(motifs))

  } else if (is.null(ArchProj) & !is.null(peaks) & class(peaks)=="GRanges") {
    if(length(peaks)==0) stop("No peaks provided.")
    message ("annotating peaks with motifs")
    BS.genome <- switch(genome,
                        hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                        hg19 = "BSgenome.Hsapiens.UCSC.hg19",
                        mm10 = "BSgenome.Mmusculus.UCSC.mm10")

    peaks.pruned <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
    if(length(peaks.pruned)==0) {
      warning("No peaks in standard chromosomes. NAs returned.")
      regulon[,field_name] <- NA
      return(regulon)
    }
    # peaks shoud by unique otherwise some idxAATAC values will be missing
    # in peaks.idx object because match function always returns the first index
    if(any(duplicated(peaks.pruned))) stop("Duplicated peaks provided.")

    # store original peak indices to match them to regulon idxATAC
    peaks.idx <- GenomicRanges::match(peaks.pruned, peaks)
    peaks.pruned <- peaks.pruned[peaks.idx %in% regulon$idxATAC]
    peaks.idx <- peaks.idx[peaks.idx %in% regulon$idxATAC]
    motifs <- epiregulon:::annotateMotif(species, peaks.pruned, BS.genome, pwms, ...)
    motifs <- assay(motifs,"motifMatches")

    # Convert motifs to gene names
    motif_names <- unlist(lapply(strsplit(colnames(motifs), split="_|\\."), "[", 3))

    colnames(motifs) <- epiregulon:::matchNames(motif_names, regulon)

  } else {

    stop("specify either an ArchR project path OR supply a GenomicRanges object for peaks")
  }



  # Remove motifs not found in regulon
  motifs <- motifs[, colnames(motifs) %in% unique(regulon$tf), drop=FALSE]

  # Add motif information
  regulon[,field_name] <- NA

  tfs_with_motif <- intersect(colnames(motifs), unique(regulon$tf))

  for (tf in tfs_with_motif){
    regulon[which(regulon$tf ==tf), field_name] <- motifs[match(regulon$idxATAC[which(regulon$tf ==tf)], peaks.idx),tf]
  }

  regulon[,field_name] <- as.numeric(regulon[,field_name])

  regulon


}



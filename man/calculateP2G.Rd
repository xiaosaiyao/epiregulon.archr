% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateP2G.R
\name{calculateP2G}
\alias{calculateP2G}
\title{Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes}
\usage{
calculateP2G(
  peakMatrix = NULL,
  expMatrix = NULL,
  reducedDim = NULL,
  ArchR_path = NULL,
  useDim = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  cor_cutoff = 0.5,
  cellNum = 100,
  exp_assay = "logcounts",
  peak_assay = "counts",
  gene_symbol = "name",
  clusters = NULL,
  cor_method = c("pearson", "kendall", "spearman"),
  ...
)
}
\arguments{
\item{peakMatrix}{A SingleCellExperiment object containing counts of chromatin accessibility at each peak region or genomic bin from scATAC-seq.
\code{rowRanges} should contain genomic positions of the peaks in the form of \code{GRanges}.}

\item{expMatrix}{A SingleCellExperiment object containing gene expression counts from scRNA-seq. \code{rowRanges} should contain genomic positions of
the genes in the form of \code{GRanges}. \code{rowData} should contain a column of gene symbols with column name matching the \code{gene_symbol} argument.}

\item{reducedDim}{A matrix of dimension reduced values}

\item{ArchR_path}{String specifying the path to an ArchR project if ArchR's implementation of addPeak2GeneLinks is desired}

\item{useDim}{String specifying the dimensional reduction representation in the ArchR project to use or the name of the reduced dimension matrix supplied by the user}

\item{useMatrix}{String specifying which the name of the gene expression matrix in the ArchR project to use.
It is often called the "GeneExpressionMatrix" for multiome and "GeneIntegrationMatrix" for unpaired data in ArchR project.}

\item{cor_cutoff}{A numeric scalar to specify the correlation cutoff between ATAC-seq peaks and RNA-seq genes to assign peak to gene links.
Default correlation cutoff is 0.5.}

\item{cellNum}{An integer to specify the number of cells to include in each K-means cluster. Default is 200 cells.}

\item{exp_assay}{String indicating the name of the assay in expMatrix for gene expression}

\item{peak_assay}{String indicating the name of the assay in peakMatrix for chromatin accessibility}

\item{gene_symbol}{String indicating the column name in the rowData of expMatrix that corresponds to gene symbol}

\item{clusters}{A vector corresponding to the cluster labels for calculation of correlations within each cluster. If left NULL, correlation is calculated across
all clusters. See details for the use of clusters}

\item{cor_method}{String indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman".}

\item{...}{other parameters to pass to \link[ArchR:addPeak2GeneLinks]{ArchR::addPeak2GeneLinks} package or to \link[epiregulon:calculateP2G]{epiregulon::calculateP2G}.}

\item{maxDist}{An integer to specify the base pair extension from transcription start start for overlap with peak regions}
}
\value{
A DataFrame of Peak to Gene correlation
}
\description{
Establish peak to gene links based on correlations between ATAC-seq peaks and RNA-seq genes
}
\details{
Cluster information is sometimes helpful to avoid the \href{https://en.wikipedia.org/wiki/Simpson\%27s_paradox}{Simpsons's paradox} in which baseline differences
between cell lines or cell types can create artificial or even inverse correlations between peak accessibility and gene expression. If Cluster information is provided,
correlation is performed within cell aggregates of each cluster.
}
\examples{
# create a mock singleCellExperiment object for gene expression matrix
set.seed(1000)
gene_sce <- scuttle::mockSCE()
gene_sce <- scuttle::logNormCounts(gene_sce)
gene_gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr3","chr4"), nrow(gene_sce)/4),
                   ranges = IRanges(start = seq(from = 1, length.out=nrow(gene_sce), by = 1000),
                   width = 100))
rownames(gene_sce) <- rownames(gene_sce)
gene_gr$name <- rownames(gene_sce)
rowRanges(gene_sce) <- gene_gr

# create a mock singleCellExperiment object for peak matrix
peak_gr <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = seq(from = 1, to = 10000, by = 1000), width = 100))
peak_counts <- matrix(sample(x = 0:4, size = ncol(gene_sce)*length(peak_gr), replace = TRUE),
                      nrow = length(peak_gr), ncol=ncol(gene_sce))
peak_sce <- SingleCellExperiment(list(counts = peak_counts), colData = colData(gene_sce))
rowRanges(peak_sce) <- peak_gr
rownames(peak_sce) <- paste0("peak",1:10)
# create a mock reducedDim matrix
reducedDim_mat <- matrix(runif(ncol(gene_sce)*50, min = 0, max = 1), nrow = ncol(gene_sce), 50)
p2g <- calculateP2G(peakMatrix = peak_sce, expMatrix = gene_sce, reducedDim = reducedDim_mat,
                    cellNum = 20, clusters = gene_sce$Treatment)
}
\author{
Xiaosai Yao, Shang-yang Chen
}

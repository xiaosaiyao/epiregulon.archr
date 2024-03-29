% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archr2MAE.R
\name{ArchRMatrix2SCE}
\alias{ArchRMatrix2SCE}
\title{Convert ArchR Gene Expression Matrix to SingleCellExperiment}
\usage{
ArchRMatrix2SCE(
  rse,
  rename = "counts",
  transform = FALSE,
  transform_method = c("log", "normalize_log"),
  log_name = "logcounts"
)
}
\arguments{
\item{rse}{a RangedSummarizedExperiment object obtained from \code{ArchR::getMatrixFromProject}}

\item{rename}{String indicating the name of the assay to rename}

\item{transform}{A logical indicating whether counts should be transformed or normalized.
See @details}

\item{transform_method}{String indicating counts should be log2 transformed or
library size normalized then log2 transformed.See @details}

\item{log_name}{String indicating the assay name for the log transformed values}
}
\value{
A SingleCellExperiment object
}
\description{
Convert ArchR Gene Expression Matrix to SingleCellExperiment
}
\details{
When gene expression matrices (either GeneIntegrationMatrix or GeneExpressionMatrix)
are retrieved from ArchR project, they are returned as \code{RangedSummarizedExperiment}. This function
converts the class of the object from \code{RangedSummarizedExperiment} to \code{SingleCellExperiment}.
During the conversion, the genomic location of the genes are transferred from \code{rowData} in the
\code{RangedSummarizedExperiment} to the \code{rowRanges} of the \code{SingleCellExperiment}.
The counts from ArchR project are library size normalized counts. This function allows the user to add
an extra assay to store log transformed counts.
}

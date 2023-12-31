% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probe2gene.R
\name{probe2gene}
\alias{probe2gene}
\alias{probe.2.gene.eset}
\title{Transformation of probe level expression to gene level expression}
\usage{
probe2gene(
  probeSE,
  chip = NA,
  from = "PROBEID",
  to = "ENTREZID",
  multi.to = "first",
  multi.from = "mean"
)
}
\arguments{
\item{probeSE}{Probe expression data.  An object of class
\code{\linkS4class{SummarizedExperiment}}.  Make sure that the
\code{\link{metadata}} contains an element named \code{annotation} that
provides the corresponding ID of a recognized platform such as
\code{hgu95av2} (Affymetrix Human Genome U95 chip).  This requires that a
corresponding \code{.db} package exists (see
\url{http://www.bioconductor.org/packages/release/BiocViews.html#___ChipName}
for available chips/packages) and that you have it installed.
Alternatively, the mapping from probe to gene can also be defined in the
\code{\link{rowData}} slot via two columns named (i) \code{PROBEID} for the
platform-specific probe ID, and (ii) \code{ENTREZID} for the corresponding
NCBI Entrez Gene ID.}

\item{chip}{Character. The ID of a recognized microarray platform. 
Only required if not provided in the \code{\link{metadata}} of \code{probeSE}
via an element named \code{annotation}.}

\item{from}{Character. ID type from which should be mapped. Corresponds to the
ID type of the names of argument \code{se}, with the default \code{PROBEID}
being appropriate if the mapping is based on Bioconductor annotation packages. 
Note that \code{from} is ignored if \code{to} is a \code{\link{rowData}} column 
of \code{probeSE}.}

\item{to}{Character. Gene ID type to which should be mapped. Corresponds to 
the gene ID type the rownames of argument \code{probeSE} should be updated with.
Note that this can also be the name of a column in the \code{\link{rowData}} 
slot of \code{probeSE} to specify user-defined mappings in which conflicts 
have been manually resolved. Defaults to \code{ENTREZID}.}

\item{multi.to}{How to resolve 1:many mappings, i.e. multiple gene IDs for a
single probe ID? This is passed on to the \code{multiVals} argument of
\code{\link{mapIds}} and can thus take several pre-defined values, but also
the form of a user-defined function. However, note that this requires that a
single gene ID is returned for each probe ID. Default is \code{"first"}, 
which accordingly returns the first gene ID mapped onto the respective probe ID.}

\item{multi.from}{How to resolve many:1 mappings, i.e. multiple probe IDs 
mapping to the same gene ID? Pre-defined options include:
\itemize{ \item 'mean' (Default): updates the respective gene expression with
 the average over the expression of all probes mapping to that gene, 
\item 'first': returns the first probe ID for each gene ID with 
multiple probe IDs,
\item 'minp' selects the probe ID with minimum p-value (according to the
\code{\link{rowData}} column \code{PVAL} of \code{probeSE}),
\item 'maxfc' selects the probe ID with maximum absolute log2 fold change 
(according to the \code{\link{rowData}} column \code{FC} of \code{probeSE}).}}
}
\value{
A \code{\linkS4class{SummarizedExperiment}} on gene level.
}
\description{
Transforms expression data on probe level to gene level expression by
summarizing all probes that are annotated to a particular gene.
}
\examples{

    # (1) reading the expression data from file
    exprs.file <- system.file("extdata/exprs.tab", package="EnrichmentBrowser")
    cdat.file <- system.file("extdata/colData.tab", package="EnrichmentBrowser")
    rdat.file <- system.file("extdata/rowData.tab", package="EnrichmentBrowser")
    probeSE <- readSE(exprs.file, cdat.file, rdat.file)
    geneSE <- probe2gene(probeSE) 

}
\seealso{
\code{\link{readSE}} for reading expression data from file,
\code{\link{deAna}} for differential expression analysis.
}
\author{
Ludwig Geistlinger
}

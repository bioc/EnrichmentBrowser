% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapIds.R
\name{idMap}
\alias{idMap}
\alias{map.ids}
\alias{idTypes}
\title{Mapping between gene ID types}
\usage{
idMap(
  obj,
  org = NA,
  from = "ENSEMBL",
  to = "ENTREZID",
  multi.to = "first",
  multi.from = "first"
)

idTypes(org)
}
\arguments{
\item{obj}{The object for which gene IDs should be mapped. Supported options
 include \itemize{
\item Gene expression dataset. An object of class 
\code{\linkS4class{SummarizedExperiment}}.
Expects the names to be of gene ID type given in argument \code{from}.
\item Gene sets. Either a list of gene sets (character vectors of gene
IDs) or a \code{\linkS4class{GeneSetCollection}} storing all gene sets.
\item Gene regulatory network. A 3-column character matrix;
1st col = IDs of regulating genes; 2nd col = IDs of regulated genes; 
3rd col = regulation effect; Use '+' and '-' for activation / inhibition.}}

\item{org}{Character. Organism in KEGG three letter code, e.g. \sQuote{hsa}
for \sQuote{Homo sapiens}.  See references.}

\item{from}{Character. Gene ID type from which should be mapped.  Corresponds
to the gene ID type of argument \code{obj}. Defaults to \code{ENSEMBL}.}

\item{to}{Character. Gene ID type to which should be mapped. Corresponds to 
the gene ID type the argument \code{obj} should be updated with.
If \code{obj} is an expression dataset of class 
\code{\linkS4class{SummarizedExperiment}}, \code{to} can also be the name of 
a column in the \code{\link{rowData}} 
slot to specify user-defined mappings in which conflicts have been 
manually resolved. Defaults to \code{ENTREZID}.}

\item{multi.to}{How to resolve 1:many mappings, i.e. multiple to.IDs for a 
single from.ID? This is passed on to the \code{multiVals} argument of 
\code{\link{mapIds}} and can thus take several pre-defined values, but also
the form of a user-defined function. However, note that this requires that a 
single to.ID is returned for each from.ID. Default is \code{"first"},
which accordingly returns the first to.ID mapped onto the respective from.ID.}

\item{multi.from}{How to resolve many:1 mappings, i.e. multiple from.IDs 
mapping to the same to.ID? Only applicable if \code{obj} is an expression 
dataset of class \code{\linkS4class{SummarizedExperiment}}. 
Pre-defined options include: 
\itemize{ \item 'first' (Default): returns the first from.ID for each to.ID
with multiple from.IDs, 
\item 'minp': selects the from.ID with minimum p-value (according to the 
\code{\link{rowData}} column \code{PVAL} of \code{obj}),
\item 'maxfc': selects the from.ID with maximum absolute log2 fold change 
(according to the \code{\link{rowData}} column \code{FC} of \code{obj}).}
Note that a user-defined function can also be supplied for custom behaviors.
This will be applied for each case where there are multiple from.IDs for a 
single to.ID, and accordingly takes the arguments \code{ids} and \code{obj}. 
The argument \code{ids} corresponds to the multiple from.IDs from which a 
single ID should be chosen, e.g. via information available in argument 
\code{obj}. See examples for a case where ids are selected based on a 
user-defined \code{\link{rowData}} column.}
}
\value{
idTypes: character vector listing the available gene ID types for
the mapping;

idMap: An object of the same class as the input argument \code{obj}, i.e.  
a \code{\linkS4class{SummarizedExperiment}} if provided an expression dataset,
a list of character vectors or a \code{\linkS4class{GeneSetCollection}} if 
provided gene sets, and a character matrix if provided a gene regulatory network.
}
\description{
Functionality to map between common gene ID types such as ENSEMBL and ENTREZ
for gene expression datasets, gene sets, and gene regulatory networks.
}
\details{
The function 'idTypes' lists the valid values which the arguments 'from'
and 'to' can take. This corresponds to the names of the available gene ID
types for the mapping.
}
\examples{

    # (1) ID mapping for gene expression datasets 
    # create an expression dataset with 3 genes and 3 samples
    se <- makeExampleData("SE", nfeat = 3, nsmpl = 3)
    names(se) <- paste0("ENSG00000000", c("003", "005", "419"))
    idMap(se, org = "hsa")

    # user-defined mapping
    rowData(se)$MYID <- c("g1", "g1", "g2")
    idMap(se, to = "MYID")    

    # data-driven resolving of many:1 mappings
    
    ## e.g. select from.ID with lowest p-value
    pcol <- configEBrowser("PVAL.COL")
    rowData(se)[[pcol]] <- c(0.001, 0.32, 0.15)
    idMap(se, to = "MYID", multi.from = "minp") 
   
    ## ... or using a customized function
    maxScore <- function(ids, se)
    {
         scores <- rowData(se)[ids, "SCORE"]
         ind <- which.max(scores)
         return(ids[ind])
    }
    rowData(se)$SCORE <- c(125.7, 33.4, 58.6)
    idMap(se, to = "MYID", multi.from = maxScore) 
           
    # (2) ID mapping for gene sets 
    # create two gene sets containing 3 genes each 
    s2 <- paste0("ENSG00000", c("012048", "139618", "141510"))
    gs <- list(s1 = names(se), s2 = s2)
    idMap(gs, org = "hsa", from = "ENSEMBL", to = "SYMBOL")    

    # (3) ID mapping for gene regulatory networks
    grn <- cbind(FROM = gs$s1, TO = gs$s2, TYPE = rep("+", 3))
    idMap(grn, org = "hsa", from = "ENSEMBL", to = "ENTREZID")  

}
\references{
KEGG Organism code
\url{http://www.genome.jp/kegg/catalog/org_list.html}
}
\seealso{
\code{\linkS4class{SummarizedExperiment}}, \code{\link{mapIds}},
\code{\link{keytypes}}
}
\author{
Ludwig Geistlinger
}

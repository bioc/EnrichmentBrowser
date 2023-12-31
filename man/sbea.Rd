% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbea.R
\name{sbeaMethods}
\alias{sbeaMethods}
\alias{sbea}
\alias{ora}
\alias{gsea}
\title{Set-based enrichment analysis (SBEA)}
\usage{
sbeaMethods()

sbea(
  method = EnrichmentBrowser::sbeaMethods(),
  se,
  gs,
  alpha = 0.05,
  perm = 1000,
  padj.method = "none",
  out.file = NULL,
  browse = FALSE,
  assay = "auto",
  ...
)
}
\arguments{
\item{method}{Set-based enrichment analysis method.  Currently, the
following set-based enrichment analysis methods are supported: \sQuote{ora},
\sQuote{safe}, \sQuote{gsea}, \sQuote{padog}, \sQuote{roast},
\sQuote{camera}, \sQuote{gsa}, \sQuote{gsva}, \sQuote{globaltest},
\sQuote{samgs}, \sQuote{ebm}, and \sQuote{mgsa}.  For basic ora also set
'perm=0'. Default is \sQuote{ora}.  This can also be a
user-defined function implementing a set-based enrichment method. See Details.}

\item{se}{Expression dataset.  An object of class
\code{\linkS4class{SummarizedExperiment}}.  Mandatory minimal annotations:
\itemize{ \item colData column storing binary group assignment (named
"GROUP") \item rowData column storing (log2) fold changes of differential
expression between sample groups (named "FC") \item rowData column storing
adjusted (corrected for multiple testing) p-values of differential
expression between sample groups (named "ADJ.PVAL") } Additional optional
annotations: \itemize{ \item colData column defining paired samples or
sample blocks (named "BLOCK") \item metadata slot named "annotation" giving
the organism under investigation in KEGG three letter code (e.g. "hsa" for
Homo sapiens) \item metadata slot named "dataType" indicating the expression
data type ("ma" for microarray, "rseq" for RNA-seq) }}

\item{gs}{Gene sets.  Either a list of gene sets (character vectors of gene
IDs) or a text file in GMT format storing all gene sets under investigation.}

\item{alpha}{Statistical significance level. Defaults to 0.05.}

\item{perm}{Number of permutations of the sample group assignments.
Defaults to 1000. For basic ora set 'perm=0'.  Using
method="gsea" and 'perm=0' invokes the permutation approximation from the
npGSEA package.}

\item{padj.method}{Method for adjusting nominal gene set p-values to
multiple testing.  For available methods see the man page of the stats
function \code{\link{p.adjust}}.  Defaults to'none', i.e. leaves the nominal
gene set p-values unadjusted.}

\item{out.file}{Optional output file the gene set ranking will be written
to.}

\item{browse}{Logical. Should results be displayed in the browser for
interactive exploration? Defaults to FALSE.}

\item{assay}{Character. The name of the assay for enrichment 
analysis if \code{se} is a \code{\linkS4class{SummarizedExperiment}} with 
*multiple assays*. Defaults to \code{"auto"}, which automatically determines
the appropriate assay based on data type provided and enrichment method selected. 
See details.}

\item{...}{Additional arguments passed to individual sbea methods.  This
includes currently for ORA and MGSA: \itemize{ \item beta: Log2 fold change
significance level. Defaults to 1 (2-fold).  \item sig.stat: decides which
statistic is used for determining significant DE genes.  Options are:
\itemize{ \item 'p' (Default): genes with adjusted p-value below alpha.  
\item 'fc': genes with abs(log2(fold change)) above beta \item '&': p & fc 
(logical AND) \item '|': p | fc (logical OR) \item 'xxp': top xx \% of genes 
sorted by adjusted p-value \item 'xxfc' top xx \% of genes sorted by absolute
log2 fold change.} }}
}
\value{
sbeaMethods: a character vector of currently supported methods;

sbea: if(is.null(out.file)): an enrichment analysis result object that can
be detailedly explored by calling \code{\link{eaBrowse}} and from which a
flat gene set ranking can be extracted by calling \code{\link{gsRanking}}.
If 'out.file' is given, the ranking is written to the specified file.
}
\description{
This is the main function for the enrichment analysis of gene sets.  It
implements and wraps existing implementations of several frequently used
methods and allows a flexible inspection of resulting gene set rankings.
}
\details{
'ora': overrepresentation analysis, simple and frequently used test based on
the hypergeometric distribution (see Goeman and Buhlmann, 2007, for a
critical review).

'safe': significance analysis of function and expression, generalization of
ORA, includes other test statistics, e.g. Wilcoxon's rank sum, and allows to
estimate the significance of gene sets by sample permutation; implemented in
the safe package (Barry et al., 2005).

'gsea': gene set enrichment analysis, frequently used and widely accepted,
uses a Kolmogorov-Smirnov statistic to test whether the ranks of the
p-values of genes in a gene set resemble a uniform distribution (Subramanian
et al., 2005).

'padog': pathway analysis with down-weighting of overlapping genes,
incorporates gene weights to favor genes appearing in few pathways versus
genes that appear in many pathways; implemented in the PADOG package.

'roast': rotation gene set test, uses rotation instead of permutation for
assessment of gene set significance; implemented in the limma and edgeR
packages for microarray and RNA-seq data, respectively.

'camera': correlation adjusted mean rank gene set test, accounts for
inter-gene correlations as implemented in the limma and edgeR packages for
microarray and RNA-seq data, respectively.

'gsa': gene set analysis, differs from GSEA by using the maxmean statistic,
i.e. the mean of the positive or negative part of gene scores in the gene
set; implemented in the GSA package.

'gsva': gene set variation analysis, transforms the data from a gene by
sample matrix to a gene set by sample matrix, thereby allowing the
evaluation of gene set enrichment for each sample; implemented in the GSVA
package.

'globaltest': global testing of groups of genes, general test of groups of
genes for association with a response variable; implemented in the
globaltest package.

'samgs': significance analysis of microarrays on gene sets, extends the SAM
method for single genes to gene set analysis (Dinu et al., 2007).

'ebm': empirical Brown's method, combines p-values of genes in a gene set
using Brown's method to combine p-values from dependent tests; implemented
in the EmpiricalBrownsMethod package.

'mgsa': model-based gene set analysis, Bayesian modeling approach taking set
overlap into account by working on all sets simultaneously, thereby reducing
the number of redundant sets; implemented in the mgsa package.

It is also possible to use additional set-based enrichment methods.  This
requires to implement a function that takes 'se' and 'gs'
as arguments and returns a numeric vector 'ps' storing the resulting p-value
for each gene set in 'gs'. This vector must be named accordingly (i.e.
names(ps) == names(gs)). See examples.

Using a \code{\linkS4class{SummarizedExperiment}} with *multiple assays*:

For the typical use case within the EnrichmentBrowser workflow this will
be a \code{\linkS4class{SummarizedExperiment}} with two assays: (i) an assay
storing the *raw* expression values, and (ii) an assay storing the *norm*alized
expression values as obtained with the \code{\link{normalize}} function. 

In this case, \code{assay = "auto"} will *auto*matically determine the assay 
based on the data type provided and the enrichment method selected. 
For usage outside of the typical workflow, the \code{assay} argument can be
used to provide the name of the assay for the enrichment analysis.
}
\examples{

    # currently supported methods
    sbeaMethods()

    # (1) expression data: 
    # simulated expression values of 100 genes
    # in two sample groups of 6 samples each
    se <- makeExampleData(what="SE")
    se <- deAna(se)

    # (2) gene sets:
    # draw 10 gene sets with 15-25 genes
    gs <- makeExampleData(what="gs", gnames=names(se))

    # (3) make 2 artificially enriched sets:
    sig.genes <- names(se)[rowData(se)$ADJ.PVAL < 0.1]
    gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
    gs[[2]] <- sample(sig.genes, length(gs[[2]]))   

    # (4) performing the enrichment analysis
    ea.res <- sbea(method="ora", se=se, gs=gs, perm=0)

    # (5) result visualization and exploration
    gsRanking(ea.res)

    # using your own tailored function as enrichment method
    dummySBEA <- function(se, gs)
    {
        sig.ps <- sample(seq(0, 0.05, length=1000), 5)
        nsig.ps <- sample(seq(0.1, 1, length=1000), length(gs)-5)
        ps <- sample(c(sig.ps, nsig.ps), length(gs))
        names(ps) <- names(gs)
        return(ps)
    }

    ea.res2 <- sbea(method=dummySBEA, se=se, gs=gs)
    gsRanking(ea.res2) 

}
\references{
Geistlinger at al. (2020) Towards a gold standard for benchmarking  
gene set enrichment analysis. Briefings in Bioinformatics.

Goeman and Buhlmann (2007) Analyzing gene expression data in
terms of gene sets: methodological issues. Bioinformatics, 23:980-7.

Subramanian et al. (2005) Gene Set Enrichment Analysis: a knowledge-based
approach for interpreting genome-wide expression profiles. PNAS, 102:15545-50.
}
\seealso{
Input: \code{\link{readSE}}, \code{\link{probe2gene}}
\code{\link{getGenesets}} to retrieve gene sets from databases such as GO 
and KEGG.

Output: \code{\link{gsRanking}} to retrieve the ranked list of gene sets.
\code{\link{eaBrowse}} for exploration of resulting gene sets.

Other: \code{\link{nbea}} to perform network-based enrichment analysis.
\code{\link{combResults}} to combine results from different methods.
}
\author{
Ludwig Geistlinger
}

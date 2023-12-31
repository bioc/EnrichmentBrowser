############################################################
#
# author: Ludwig Geistlinger
# date: 3 Feb 2011
#
# GGEA - Gene Graph Enrichment Analysis
#
# Update, 09 May 2014: Extension to network-based enrichment
#           analysis 
#
############################################################

#' @rdname nbea
#' @export
nbeaMethods <- function() 
    c("ggea", "spia", "pathnet", "degraph",
	  "ganpa", "cepa", "topologygsa", "netgsa", "neat")


#' Network-based enrichment analysis (NBEA)
#' 
#' This is the main function for network-based enrichment analysis.  It
#' implements and wraps existing implementations of several frequently used
#' methods and allows a flexible inspection of resulting gene set rankings.
#' 
#' 'ggea': gene graph enrichment analysis, scores gene sets according to
#' consistency within the given gene regulatory network, i.e. checks activating
#' regulations for positive correlation and repressing regulations for negative
#' correlation of regulator and target gene expression (Geistlinger et al.,
#' 2011). When using 'ggea' it is possible to estimate the statistical
#' significance of the consistency score of each gene set in two different
#' ways: (1) based on sample permutation as described in the original
#' publication (Geistlinger et al., 2011) or (2) using an approximation in the
#' spirit of Bioconductor's npGSEA package that is much faster.
#' 
#' 'spia': signaling pathway impact analysis, combines ORA with the probability
#' that expression changes are propagated across the pathway topology;
#' implemented in Bioconductor's SPIA package (Tarca et al., 2009).
#' 
#' 'pathnet': pathway analysis using network information, applies ORA on
#' combined evidence for the observed signal for gene nodes and the signal
#' implied by connected neighbors in the network; implemented in Bioconductor's
#' PathNet package.
#' 
#' 'degraph': differential expression testing for gene graphs, multivariate
#' testing of differences in mean incorporating underlying graph structure;
#' implemented in Bioconductor's DEGraph package.
#' 
#' 'topologygsa': topology-based gene set analysis, uses Gaussian graphical
#' models to incorporate the dependence structure among genes as implied by
#' pathway topology; implemented in CRAN's topologyGSA package.
#' 
#' 'ganpa': gene association network-based pathway analysis, incorporates
#' network-derived gene weights in the enrichment analysis; implemented in
#' CRAN's GANPA package.
#' 
#' 'cepa': centrality-based pathway enrichment, incorporates network
#' centralities as node weights mapped from differentially expressed genes in
#' pathways; implemented in CRAN's CePa package.
#' 
#' 'netgsa': network-based gene set analysis, incorporates external information
#' about interactions among genes as well as novel interactions learned from
#' data; implemented in CRAN's NetGSA package.
#' 
#' 'neat': network enrichment analysis test, compares the number of links between
#' differentially expressed genes and a gene set of interest to the number of
#' links expected under a hypergeometric null model; proposed by Signorelli et al.
#' (2016) and implemented in CRAN's neat package.
#' 
#' It is also possible to use additional network-based enrichment methods.
#' This requires to implement a function that takes 'se', 'gs', and 'grn'
#' and as arguments and returns a numeric matrix 'res.tbl' with a
#' mandatory column named 'PVAL' storing the resulting p-value for each gene
#' set in 'gs'. The rows of this matrix must be named accordingly (i.e.
#' rownames(res.tbl) == names(gs)). See examples.
#' 
#' Using a \code{\linkS4class{SummarizedExperiment}} with *multiple assays*:
#' 
#' For the typical use case within the EnrichmentBrowser workflow this will
#' be a \code{\linkS4class{SummarizedExperiment}} with two assays: (i) an assay
#' storing the *raw* expression values, and (ii) an assay storing the *norm*alized
#' expression values as obtained with the \code{\link{normalize}} function. 
#' 
#' In this case, \code{assay = "auto"} will *auto*matically determine the assay 
#' based on the data type provided and the enrichment method selected. 
#' For usage outside of the typical workflow, the \code{assay} argument can be
#' used to provide the name of the assay for the enrichment analysis.
#'
#' @aliases ggea spia
#' @param method Network-based enrichment analysis method.  Currently, the
#' following network-based enrichment analysis methods are supported:
#' \sQuote{ggea}, \sQuote{spia}, \sQuote{pathnet}, \sQuote{degraph},
#' \sQuote{topologygsa}, \sQuote{ganpa}, \sQuote{cepa}, \sQuote{netgsa},
#' and \sQuote{neat}. Default is 'ggea'.  This can also be the name of a
#' user-defined function implementing a method for network-based enrichment analysis.
#' See Details.
#' @param se Expression dataset.  An object of class
#' \code{\linkS4class{SummarizedExperiment}}.  Mandatory minimal annotations:
#' \itemize{ \item colData column storing binary group assignment (named
#' "GROUP") \item rowData column storing (log2) fold changes of differential
#' expression between sample groups (named "FC") \item rowData column storing
#' adjusted (corrected for multiple testing) p-values of differential
#' expression between sample groups (named "ADJ.PVAL") } Additional optional
#' annotations: \itemize{ \item colData column defining paired samples or
#' sample blocks (named "BLOCK") \item metadata slot named "annotation" giving
#' the organism under investigation in KEGG three letter code (e.g. "hsa" for
#' Homo sapiens) \item metadata slot named "dataType" indicating the expression
#' data type ("ma" for microarray, "rseq" for RNA-seq) }
#' @param gs Gene sets.  Either a list of gene sets (character vectors of gene
#' IDs) or a text file in GMT format storing all gene sets under investigation.
#' @param grn Gene regulatory network.  Either an absolute file path to a
#' tabular file or a character matrix with exactly *THREE* cols; 1st col = IDs
#' of regulating genes; 2nd col = corresponding regulated genes; 3rd col =
#' regulation effect; Use '+' and '-' for activation/inhibition.
#' @param prune.grn Logical.  Should the GRN be pruned? This removes
#' duplicated, self, and reversed edges. Defaults to TRUE.
#' @param alpha Statistical significance level. Defaults to 0.05.
#' @param perm Number of permutations of the expression matrix to estimate the
#' null distribution. Defaults to 1000. If using method=\sQuote{ggea}, it is
#' possible to set 'perm=0' to use a fast approximation of gene set
#' significance to avoid permutation testing. See Details.
#' @param padj.method Method for adjusting nominal gene set p-values to
#' multiple testing.  For available methods see the man page of the stats
#' function \code{\link{p.adjust}}.  Defaults to 'none', i.e. leaves the
#' nominal gene set p-values unadjusted.
#' @param out.file Optional output file the gene set ranking will be written
#' to.
#' @param browse Logical. Should results be displayed in the browser for
#' interactive exploration? Defaults to FALSE.
#' @param assay Character. The name of the assay for enrichment 
#' analysis if \code{se} is a \code{\linkS4class{SummarizedExperiment}} with 
#' *multiple assays*. Defaults to \code{"auto"}, which automatically determines
#' the appropriate assay based on data type provided and enrichment method selected. 
#' See details.   
#' @param ...  Additional arguments passed to individual nbea methods.  This
#' includes currently: \itemize{ \item beta: Log2 fold change significance
#' level. Defaults to 1 (2-fold).  } For SPIA and NEAT: \itemize{ \item
#' sig.stat: decides which statistic is used for determining significant DE
#' genes.  Options are: \itemize{ \item 'p' (Default): genes with adjusted p-value below
#' alpha.  \item 'fc': genes with abs(log2(fold change)) above beta \item '&':
#' p & fc (logical AND) \item '|': p | fc (logical OR) \item 'xxp': top xx \% of genes 
#' sorted by adjusted p-value \item 'xxfc' top xx \% of genes sorted by absolute
#' log2 fold change.} } For GGEA: \itemize{
#' \item cons.thresh: edge consistency threshold between -1 and 1.  Defaults to
#' 0.2, i.e. only edges of the GRN with consistency >= 0.2 are included in the
#' analysis. Evaluation on real datasets has shown that this works best to
#' distinguish relevant gene sets.  Use consistency of -1 to include all edges.
#' \item gs.edges: decides which edges of the grn are considered for a gene set
#' under investigation. Should be one out of c('&', '|'), denoting logical AND
#' and OR. respectively. Accordingly, this either includes edges for which
#' regulator AND / OR target gene are members of the investigated gene set.  }
#' @return nbeaMethods: a character vector of currently supported methods;
#' 
#' nbea: if(is.null(out.file)): an enrichment analysis result object that can
#' be detailedly explored by calling \code{\link{eaBrowse}} and from which a
#' flat gene set ranking can be extracted by calling \code{\link{gsRanking}}.
#' If 'out.file' is given, the ranking is written to the specified file.
#' @author Ludwig Geistlinger
#' @seealso Input: \code{\link{readSE}}, \code{\link{probe2gene}},
#' \code{\link{getGenesets}} to retrieve gene set definitions from databases 
#' such as GO and KEGG.
#' \code{\link{compileGRN}} to construct a GRN from pathway databases.
#' 
#' Output: \code{\link{gsRanking}} to rank the list of gene sets.
#' \code{\link{eaBrowse}} for exploration of resulting gene sets.
#' 
#' Other: \code{\link{sbea}} to perform set-based enrichment analysis.
#' \code{\link{combResults}} to combine results from different methods.
#' @references Geistlinger et al. (2011) From sets to graphs: towards a
#' realistic enrichment analysis of transcriptomic systems. Bioinformatics,
#' 27(13), i366-73.
#'
#' Tarca et al. (2009) A novel signaling pathway impact analysis. 
#' Bioinformatics, 25(1):75-82.
#' 
#' Signorelli et al. (2016) NEAT: an efficient network enrichment analysis test.
#' BMC Bioinformatics, 17:352.
#'
#' @examples
#' 
#'     # currently supported methods
#'     nbeaMethods()
#' 
#'     # (1) expression data: 
#'     # simulated expression values of 100 genes
#'     # in two sample groups of 6 samples each
#'     se <- makeExampleData(what="SE")
#'     se <- deAna(se)
#' 
#'     # (2) gene sets:
#'     # draw 10 gene sets with 15-25 genes
#'     gs <- makeExampleData(what="gs", gnames=names(se))
#' 
#'     # (3) make 2 artificially enriched sets:
#'     sig.genes <- names(se)[rowData(se)$ADJ.PVAL < 0.1]
#'     gs[[1]] <- sample(sig.genes, length(gs[[1]])) 
#'     gs[[2]] <- sample(sig.genes, length(gs[[2]]))   
#'    
#'     # (4) gene regulatory network 
#'     grn <- makeExampleData(what="grn", nodes=names(se))
#'     
#'     # (5) performing the enrichment analysis
#'     ea.res <- nbea(method="ggea", se=se, gs=gs, grn=grn)
#' 
#'     # (6) result visualization and exploration
#'     gsRanking(ea.res, signif.only=FALSE)
#' 
#'     # using your own function as enrichment method
#'     dummyNBEA <- function(se, gs, grn)
#'     {
#'         sig.ps <- sample(seq(0,0.05, length=1000),5)
#'         insig.ps <- sample(seq(0.1,1, length=1000), length(gs)-5)
#'         ps <- sample(c(sig.ps, insig.ps), length(gs))
#'         score <- sample(1:100, length(gs), replace=TRUE)
#'         res.tbl <- cbind(score, ps)
#'         colnames(res.tbl) <- c("SCORE", "PVAL")
#'         rownames(res.tbl) <- names(gs)
#'         return(res.tbl[order(ps),])
#'     }
#' 
#'     ea.res2 <- nbea(method=dummyNBEA, se=se, gs=gs, grn=grn)
#'     gsRanking(ea.res2) 
#' 
#' @export nbea
nbea <- function(
    method=EnrichmentBrowser::nbeaMethods(), 
    se, 
    gs, 
    grn,
    prune.grn=TRUE,
    alpha=0.05, 
    perm=1000, 
    padj.method="none",
    out.file=NULL,
    browse=FALSE,
    assay = "auto", 
    ...)
{
    # get configuration
    GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
    PVAL.COL <- configEBrowser("PVAL.COL")
    ADJP.COL <-  configEBrowser("ADJP.COL")

    # TODO: disentangle DE and EA analysis
    se <- .preprocSE(se)
    se <- .setAssay(method, se, perm, assay)
    
    # getting gene sets & grn
    if(is(gs, "GeneSetCollection")) gs <- GSEABase::geneIds(gs)
    if(!is.list(gs)) gs <- getGenesets(gs)
    if(!is.matrix(grn)) grn <- .readGRN(grn)

    # prune grn
    if(prune.grn) grn <- .pruneGRN(grn)

    # restrict to relevant genes 
    # in the intersection of se, gs, and grn
    gs.genes <- unique(unlist(gs))
    grn.genes <- unique(c(grn[,1], grn[,2]))
    se.genes <- rownames(se)
    rel.genes <- intersect(intersect(gs.genes, grn.genes), se.genes)
    if(!length(rel.genes)) 
        stop(paste("Expression dataset (se), gene sets (gs), and", 
                    "gene regulatory network (grn) have no gene IDs in common"))

    se <- se[rel.genes,]
    gs <- lapply(gs, function(s) s[s%in% rel.genes])
    lens <- lengths(gs)
    gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
    grn <- grn[grn[,1] %in% rel.genes & grn[,2] %in% rel.genes,] 
    
    # execute ea
    if(is.character(method))
    {
        method <- match.arg(method)
        if(method == "spia") res.tbl <- .spia(se, gs, grn, alpha, perm, ...)
        else if(method == "pathnet") res.tbl <- .pathnet(se, gs, grn, alpha, perm)
        else if(method == "netgsa") res.tbl <- .netgsa(se, gs, grn)
        else if(method == "ganpa") res.tbl <- .ganpa(se, gs, grn, perm)
        else if(method == "cepa") res.tbl <- .cepa(se, gs, grn)
        else if(method == "degraph") res.tbl <- .degraph(se, gs, grn)
        else if(method == "topologygsa") res.tbl <- .topogsa(se, gs, grn, alpha, perm, ...)
        else if(method == "neat") res.tbl <- .neat(se, gs, grn, alpha, ...)
        else res.tbl <- .ggea(se, gs, grn, alpha, perm=perm, ...)      
    }
    else if(is.function(method)) 
        res.tbl <- method(se=se, gs=gs, grn=grn, ...)
    else stop(paste(method, "is not a valid method for nbea"))

    res.tbl <- .formatEAResult(res.tbl, padj.method, out.file)

    pcol <- ifelse(padj.method == "none", PVAL.COL, ADJP.COL)
    res <- list(
            method=method, res.tbl=res.tbl, 
            nr.sigs=sum(res.tbl[,pcol] < alpha),
            se=se, gs=gs, alpha=alpha)

    if(browse) eaBrowse(res)
    return(res)
}

#
# general helpers
#
.pruneGRN <- function(grn)
{
    # remove duplicates
    grn <- grn[!duplicated(grn[,1:2]),]

    # rm self edges
    grn <- grn[grn[,1] != grn[,2],]
   
    # rm rev edges 
    genes <- unique(as.vector(grn[,1:2]))
    ggrid <- seq_along(genes)
    names(ggrid) <- genes
    igrn <- cbind(ggrid[grn[,1]], ggrid[grn[,2]])

    n <- nrow(grn)
    grid <- seq_len(n-1)
    ind <- vapply(grid,
        function(i)
        {
            x <- igrn[i,2:1]
            j <- i + 1
            cigrn <- igrn[j:n,,drop=FALSE]
            cigrn <- cigrn[cigrn[,1] == x[1], , drop=FALSE]
            is.rev <- any( cigrn[,2] == x[2] )
            return(is.rev)
        }, logical(1))
    ind <- c(ind, FALSE)
    grn <- grn[!ind,]
}
    
#
# 1 SPIA
#
.spia <- function(se, gs, grn, 
    alpha=0.05, perm=1000, beta=1, sig.stat=c("p", "fc", "|", "&")) 
{
    FC.COL <- configEBrowser("FC.COL")
    ADJP.COL <- configEBrowser("ADJP.COL")
    PVAL.COL <- configEBrowser("PVAL.COL")

    de.genes <- .isSig(rowData(se), alpha, beta, sig.stat)
    de <- rowData(se)[de.genes, FC.COL]
    names(de) <- rownames(se)[de.genes]
    all <- rownames(se)

    is.kegg <- .detectGSType(names(gs)[1]) == "KEGG"
    organism <- ""
    data.dir <- NULL
    if(is.kegg)
    { 
        id <- unlist(strsplit(names(gs)[1], "_"))
        organism <- sub("[0-9]+$", "", id[1])
    }
    else
    {     
        message("making SPIA data ...")
        path.info <- .makeSPIAData(gs, grn)
        data.dir <- system.file("extdata", package="SPIA")
        save(path.info, file=file.path(data.dir, "SPIA.RData"))
        data.dir <- paste0(data.dir, "/")
    }

    res <- SPIA::spia(de=de, all=all, organism=organism, data.dir=data.dir, nB=perm)
    
    if(is.kegg)
    {
        spl <- strsplit(names(gs), "_")
        ids <- vapply(spl, function(n) n[1], character(1))
        ids <- sub(organism, "", ids)   
        res <- res[res$ID %in% ids, ]
        names(spl) <- ids
        .acc <- function(s) paste(s[2:length(s)], collapse = " ")
        res$Name <- vapply(spl[res$ID], .acc, character(1))
    }

    res[,"Name"] <- gsub(" ", "_", res[,"Name"])
    rownames(res) <- paste(paste0(organism, res[,"ID"]), res[,"Name"], sep="_")
    res <- res[, c("pSize", "NDE", "tA", "Status", "pG")]
    colnames(res) <- c("SIZE", "NDE", "T.ACT", "STATUS", PVAL.COL)
    res[,"STATUS"] <- ifelse(res[,"STATUS"] == "Activated", 1, -1)
    as.matrix(res)
}

# spia helper: create pathway data from gs and grn
.makeSPIAData <- function(gs, grn)
{
    rel <- c("activation", "compound", "binding/association", 
            "expression", "inhibition", "activation_phosphorylation", 
            "phosphorylation", "inhibition_phosphorylation", 
            "inhibition_dephosphorylation", "dissociation", "dephosphorylation", 
            "activation_dephosphorylation", "state change", "activation_indirect effect", 
            "inhibition_ubiquination", "ubiquination", "expression_indirect effect", 
            "inhibition_indirect effect", "repression", "dissociation_phosphorylation", 
            "indirect effect_phosphorylation", "activation_binding/association", 
            "indirect effect", "activation_compound", "activation_ubiquination")

    spia.data <- sapply(names(gs), 
        function(s)
        {   
            x <- gs[[s]]
            len <- length(x) 
            nam <- list(x, x)
            m <- matrix(0, nrow=len, ncol=len, dimnames=nam)
            sgrn <- .queryGRN(gs=x, grn=grn, index=FALSE)
            if(nrow(sgrn) < configEBrowser("GS.MIN.SIZE")) return(NULL)
            act.grn <- sgrn[sgrn[,3] == "+",,drop=FALSE]
            actm2 <- m
            if(nrow(act.grn))
            {
                if(nrow(act.grn) > 1) actm <- .grn2adjm(act.grn)
                else actm <- matrix(1, nrow=1, ncol=1, dimnames=list(act.grn[1,1], act.grn[1,2]))
                actm2[rownames(actm), colnames(actm)] <- actm
            }
            inh.grn <- sgrn[sgrn[,3] == "-",,drop=FALSE]
            inhm2 <- m
            if(nrow(inh.grn))
            {
                if(nrow(inh.grn) > 1) inhm <- .grn2adjm(inh.grn)
                else inhm <-  matrix(1, nrow=1, ncol=1, dimnames=list(inh.grn[1,1], inh.grn[1,2]))
                inhm2[rownames(inhm), colnames(inhm)] <- inhm
            }
            l <- lapply(rel, 
                function(r)
                {
                    if(r == "activation") return(actm2)
                    else if(r == "inhibition") return(inhm2)
                    else return(m)                    
                } 
            )
            names(l) <- rel    
            l$nodes <- x
            l$title <- s
            l$NumberOfReactions <- 0
            return(l)
        }
    )
    spia.data <- spia.data[!sapply(spia.data, is.null)]
    return(spia.data)
}

#
# 2 NEA
#
.nea <- function(se, gs, grn, 
    alpha=0.05, perm=100, beta=1, sig.stat=c("p", "fc", "|", "&"))
{
    nea <- NULL
    isAvailable("neaGUI", type="software")

    #if(perm > 100) perm <- 100
    isig <- .isSig(rowData(se), alpha, beta, sig.stat)
    ags <- rownames(se)[isig]
    grn <- unique(grn[,1:2])
    gs.genes <- unique(unlist(gs))
    grn <- grn[(grn[,1] %in% gs.genes) & (grn[,2] %in% gs.genes),]
    network <- apply(grn, 1, function(x) paste(x, collapse=" "))
    message("Computing NEA permutations, this may take a few minutes ...")
    res <- nea(ags=ags, fgs=gs, network=network, nperm=perm)
    res <- res$MainResult
    res <- res[, c("Number_of_Genes", 
        "Number_of_AGS_genes", "Number_links", "Z_score", "P_value")]
    res <- res[order(res[,"Z_score"], decreasing=TRUE), ]
    colnames(res) <- sub("AGS", "de", colnames(res))
    colnames(res) <- sub("Number", "Nr", colnames(res))
    colnames(res) <- sub("_of", "", colnames(res))
    colnames(res) <- gsub("_", ".", colnames(res))
    colnames(res) <- toupper(colnames(res))
    res <- as.matrix(res)
    PVAL.COL <- configEBrowser("PVAL.COL")
    res <- res[order(res[,PVAL.COL]),]
    return(res) 
}

#
# 3 Pathnet
#
.pathnet <- function(se, gs, grn, alpha=0.05, perm=1000)
{
    PathNet <- NULL
    isAvailable("PathNet", type="software")
    
    ADJP.COL <- configEBrowser("ADJP.COL")
    PVAL.COL <- configEBrowser("PVAL.COL")

    dir.evid <- -log(rowData(se)[,ADJP.COL], base=10)
    dir.evid <- cbind(as.integer(rownames(se)), dir.evid)
    colnames(dir.evid) <- c("Gene.ID", "Obs")
    adjm <- .grn2adjm(grn, directed=FALSE)
    pwy.dat <- .extrPwyDat(gs, grn)
    
    res <- PathNet(
            #Enrichment_Analysis = TRUE, 
            #Contextual_Analysis = FALSE, 
            DirectEvidence_info = dir.evid, 
            Column_DirectEvidence = 2,
            Adjacency = adjm, 
            pathway = pwy.dat, 
            n_perm = perm, 
            threshold = alpha)#,
            #use_sig_pathways  = FALSE)

    res <- res$enrichment_results[, 
        c("Name", "No_of_Genes", "Sig_Direct", "Sig_Combi", "p_PathNet")]
    rownames(res) <- vapply(as.vector(res[,1]), 
        function(s) grep(paste0("^", unlist(strsplit(s,"_"))[1], "_"), names(gs), value=TRUE),
        character(1), USE.NAMES=FALSE)
    res <- res[-1]    
    colnames(res) <- c("NR.GENES", "NR.SIG.GENES", "NR.SIG.COMB.GENES", PVAL.COL)
    as.matrix(res)
}

# pathnet helper: extract pathway data from gs and grn
.extrPwyDat <- function(gs, grn)
{
    pwy.dat <- lapply(names(gs), 
        function(n)
        {
            genes <- gs[[n]] 
            sgrn <- .queryGRN(gs=genes, grn=grn, index=FALSE)
            if(nrow(sgrn))
                dat <- cbind(sgrn[,1:2, drop=FALSE], rep(n, nrow(sgrn)))
            else dat <- NULL
        }
    )
    ind <- vapply(pwy.dat, is.null, logical(1))
    pwy.dat <- pwy.dat[!ind]
    nr <- sum( vapply(pwy.dat, nrow, integer(1)) )
    pwy.datm <- matrix("", nrow=nr, ncol=3)
    colnames(pwy.datm) <- c("id1", "id2", "title")
    start <- 1
    for(i in seq_len(length(pwy.dat)))
    {
        end <- start + nrow(pwy.dat[[i]]) - 1
        pwy.datm[start:end,] <- pwy.dat[[i]]
        start <- end + 1
    }
    pwy.dat <- data.frame(id1=as.integer(pwy.datm[,1]), 
        id2=as.integer(pwy.datm[,2]), title=pwy.datm[,3])
    return(pwy.dat)
}

# pathnet helper: converts 3-col grn to adjacency matrix
.grn2adjm <- function(grn, directed=TRUE)
{
    nodes <- sort(unique(as.vector(grn[,1:2])))
    adjm <- sapply(nodes, 
        function(n)
        {
            tgs <- grep(n, grn[,1])
            if(length(tgs))
            {
                tgs <- grn[tgs,2]
                adjv <- as.integer(nodes %in% tgs)
            }
            else adjv <- rep(0, length(nodes))
            return(adjv) 
        }) 
    rownames(adjm) <- nodes
    adjm <- t(adjm)
  
    if(!directed)
        for(i in seq_along(nodes)) 
            for(j in seq_along(nodes)) 
                if(adjm[i,j]) adjm[j,i] <- 1

    return(adjm)
}

#
# 4 NetGSA
#
.netgsa <- function(se, gs, grn)
{
    NetGSA <- bic.netEst.undir <- netEst.undir <- prepareAdjacencyMatrix <- NULL
    isAvailable("netgsa", type="software")

    x <- assay(se)
    y <- colData(se)[,configEBrowser("GRP.COL")] + 1

    # prepare network
    out.dir <- configEBrowser("OUTDIR.DEFAULT")
    if(!file.exists(out.dir)) dir.create(out.dir, recursive=TRUE)
    out.file <- file.path(out.dir, "grn.txt")
    grn <- cbind(grn[,1:2], "undirected")
    colnames(grn) <- c("src", "dest", "direction")
    write.csv(grn, file=out.file, row.names=FALSE)
    
    adj <- prepareAdjacencyMatrix(x, y, gs, FALSE, out.file, NULL)
    grn <- adj$Adj
    cmat <- adj$B
    
    file.remove(out.file)
    ind <- intersect(rownames(x), rownames(grn))
    grn <- grn[ind, ind]
    x <- x[rownames(grn),]
    cmat <- cmat[,rownames(grn)]
    cmat <- cmat[rowSums(cmat) > 2,]
   
    message("Estimating weighted adjacency matrix for GRN") 
    message("This may take a while ...")
    .fitM <- function(i)
    {
        dat <- x[,y == i]
        lambda <- sqrt(log(nrow(dat)) / ncol(dat))
        fit.BIC <- bic.netEst.undir(dat, one = grn, lambda = lambda, eta = 0.1)
        lambda <- which.min(fit.BIC$BIC) * lambda
        fit <- netEst.undir(dat, one = grn, lambda = lambda, eta = 0.1)
        fit$Adj
    }
    A <- lapply(1:2, .fitM) 

    # execute
    message("Executing NetGSA ...")
    res <- NetGSA(A, x, y, cmat, lklMethod = 'REHE')

    res <- res$results
    rownames(res) <- res[,"pathway"]
    res <- res[,c("pSize", "pval")]
    colnames(res) <- c("NR.GENES", configEBrowser("PVAL.COL"))
    return(res)
}

#
# 5 GANPA
#
.ganpa <- function(se, gs, grn, perm=1000)
{
    GSE.Test.Main <- NULL
    isAvailable("GANPA", type="software")

    # configure
    GRP.COL <- configEBrowser("GRP.COL")
    SMPL.COL <- configEBrowser("SMPL.COL")
    OUT.DIR <- configEBrowser("OUTDIR.DEFAULT")
    GS.MIN.SIZE <- configEBrowser("GS.MIN.SIZE")
    GS.MAX.SIZE <- configEBrowser("GS.MAX.SIZE")
    PVAL.COL <- configEBrowser("PVAL.COL")
    
    if(!file.exists(OUT.DIR)) dir.create(OUT.DIR, recursive=TRUE)
    out.prefix <- file.path(OUT.DIR, "ganpa")

    # expression data
    has.scol <- SMPL.COL %in% colnames(colData(se))
    if(!has.scol) colData(se)[,SMPL.COL] <- colnames(se)
    sinfo <- colData(se)[,c(SMPL.COL, GRP.COL)]
    colnames(sinfo) <- c("sampleid", "status")
    expr.obj <- list(gExprs=assay(se), sampleinfo=sinfo)
    
    # gene regulatory network
    gnet <- .grn2gnet(grn)    

    # execute
    GSE.Test.Main(gExprs.obj=expr.obj, gsets=gs, gNET=gnet, 
        permN=perm, size.min=GS.MIN.SIZE, size.max=GS.MAX.SIZE,
        msp.correction=FALSE, output.label=out.prefix, permFDR.cutoff=1)

    # read results from output csv
    res <- read.csv(paste(out.prefix, "MeanAbs.OrigW.csv", sep=".")) 
    n <- res[,1]
    res <- as.matrix(res[,c("Size", "S", "NS", "permP")])
    colnames(res)[c(1,4)] <- c("SIZE", PVAL.COL)
    rownames(res) <- n
    return(res)
}

.grn2gnet <- function(grn)
{
    ureg <- unique(grn[,1])
    gnet <- sapply(ureg, function(r) grn[grn[,1] == r,2])
    return(gnet)
}

#
# 6 CePa
#
.cepa <- function(se, gs, grn, perm=1000)
{
    cepa.all <- set.pathway.catalogue <- sampleLabel <- NULL
    isAvailable("CePa", type="software")

    # define sample groups
    GRP.COL <- configEBrowser("GRP.COL")
    sl <- sampleLabel(colData(se)[, GRP.COL], treatment=1, control=0)

    # create pathway catalogue from gs and grn
    # (1) pathway list
    pl <- sapply(gs, function(s) as.character(.queryGRN(s, grn)))
    pl <- pl[sapply(pl, length) >= configEBrowser("GS.MIN.SIZE")]

    # (2) interaction list
    il <- data.frame(as.character(seq_len(nrow(grn))), grn[,1:2], stringsAsFactors=FALSE)
    colnames(il) <- c("interaction.id", "input", "output")

    # (3) mapping
    m <- data.frame(node.id=rownames(se), symbol=rownames(se), stringsAsFactors=FALSE)
    pwy.cat <- set.pathway.catalogue(pathList=pl, interactionList=il, mapping=m)
    
    # executing
    res <- cepa.all(mat=assay(se), label=sl, pc=pwy.cat, iter=perm)
    res.mat <- matrix(0.0, nrow=length(res), ncol=7)
    for(i in seq_along(res))
    {
        res.mat[i,1:6] <- sapply(res[[i]], function(x) x$p.value)
        res.mat[i,7] <- min(6 * min(res.mat[i,1:6]), 1) 

    }
    rownames(res.mat) <- names(res)
    n <- paste(toupper(names(res[[1]])), "PVAL" , sep=".")
    colnames(res.mat) <- c(n, configEBrowser("PVAL.COL"))
    return(res.mat)
}

#
# 7 DEGraph
#
.degraph <- function(se, gs, grn)
{    
    testOneGraph <- NULL
    isAvailable("DEGraph", type="software")

    grp <- colData(se)[,configEBrowser("GRP.COL")]
            
    options(show.error.messages=FALSE) 
    res <- sapply(names(gs),
        function(s)
        {
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            r <- try( testOneGraph(graph=gr, data=assay(se), 
                       classes=grp, useInteractionSigns=FALSE), silent=TRUE )
            if(is(r, "try-error")) return(NA) else return(r[[1]]$p.value[1])
        })
    options(show.error.messages=TRUE)
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 8 topologyGSA 
#
.topogsa <- function(se, gs, grn, alpha=0.05, perm=1000, test.mean=TRUE)
{    
    # call topologyGSA via clipper's pathQ function
    return(.clipper(se, gs, grn, alpha, perm))

    # original topologyGSA: deprecated
    # does not terminate on particular gs, eg. hsa04060_Cytokine-cytokine_receptor_interaction
    pathway.mean.test <- pathway.var.test <- NULL
    isAvailable("topologyGSA", type="software")
  
    is.DAG <- NULL
    isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     isAvailable("igraph", type="software")
 
    grp <- colData(se)[,configEBrowser("GRP.COL")]
    y1 <- t(assay(se)[, grp == 0])
    y2 <- t(assay(se)[, grp == 1])

    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            if(test.mean) r <- pathway.mean.test(y1, y2, gr, alpha, perm)
            else r <- pathway.var.test(y1, y2, gr, alpha)
            return(r$p.value)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 9 clipper
#
.clipper <- function(se, gs, grn, alpha=0.05, perm=1000)
{    
    pathQ <- NULL
    isAvailable("clipper", type="software")
  
    is.DAG <- NULL
    isAvailable("gRbase", type="software")
    
    graph_from_graphnel <- mst <- NULL
     isAvailable("igraph", type="software")
 
    grp <- colData(se)[,configEBrowser("GRP.COL")] + 1

    res <- sapply(names(gs),
        function(s)
        {
            message(s)
            gs.grn <- .queryGRN(gs[[s]], grn, index=FALSE)
            if(nrow(gs.grn) < configEBrowser("GS.MIN.SIZE")) return(NA)
            am <- .grn2adjm(gs.grn)
            gr <- graphAM(am, "directed")           
            gr <- as(gr, "graphNEL")                 
            
            if(!is.DAG(gr))
            {
                gr2 <- graph_from_graphnel(gr)
                gr2 <- mst(gr2)
                gr <- as(gr2, "graphNEL")
            }
            r <- try(pathQ(assay(se), grp, gr, perm, alpha), silent=TRUE)
            if(is(r, "try-error")) return(NA) else return(r$alphaMean)
        }
    )
    names(res) <- names(gs)
    res <- res[!is.na(res)]
    return(res)
}

#
# 10
#
.neat <- function(se, gs, grn, 
                  alpha, beta = 1, 
                  sig.stat = c("p", "fc", "|", "&"),
                  directed = TRUE)
{
    neat <- NULL
    isAvailable("neat", type = "software")

    # derive genes for alist and blist
    isig <- .isSig(rowData(se), alpha, beta, sig.stat)
    de.list <- list(de = rownames(se)[isig])
  
    # get network inputs
    grn <- grn[,1:2]
    all.nodes <- unique(as.vector(grn))

    # execute neat
    nettype <- ifelse(directed, "directed", "undirected")
    res <- neat(alist = de.list, blist = gs, network = grn, 
                nettype = nettype, nodes = all.nodes, mtc.type = "none")

    # restructure output of neat
    rownames(res) <- res$B
    res <- res[,c("nab", "expected_nab", "pvalue")]
    res <- as.matrix(res)
    colnames(res) <- c("OBS.LINKS", "EXP.LINKS", configEBrowser("PVAL.COL"))
    return(res) 
}

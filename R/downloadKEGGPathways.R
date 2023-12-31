# 
# Author: Ludwig Geistlinger
# Date: 5th August 2010
#
# download all pathways of a particular organism from the KEGG FTP
#
# UPDATE: 16 Jan 2014 - switch to KEGGREST functionality as KEGG FTP 
#           is no longer accessible
#
###############################################################################

#' Download of KEGG pathways for a particular organism
#' 
#' The function downloads all metabolic and non-metabolic pathways in KEGG XML
#' format for a specified organism.
#' 
#' 
#' @aliases download.kegg.pathways
#' @param org Organism in KEGG three letter code, e.g. \sQuote{hsa} for
#' \sQuote{homo sapiens}.
#' @param cache Logical.  Should a locally cached version used if available?
#' Defaults to \code{TRUE}.
#' @param out.dir Output directory.  If not null, pathways are written to files
#' in the specified directory.
#' @param zip Logical.  In case pathways are written to file (\sQuote{out.dir}
#' is not null): should output files be zipped?
#' @return if(is.null(out.dir)): a list of KEGGPathway objects else: none, as
#' pathways are written to file
#' @author Ludwig Geistlinger
#' @seealso \code{\link{keggList}}, \code{\link{keggGet}},
#' \code{\linkS4class{KEGGPathway}}, \code{\link{parseKGML}}
#' @examples
#' 
#'     \donttest{
#'         pwys <- downloadPathways("hsa")
#'     }
#' 
#' @export downloadPathways
downloadPathways <- function(org, cache=TRUE, out.dir=NULL, zip=FALSE)
{
    # should a cached version be used?
    pwys.name <- paste(org, "kegg", "pwys", sep=".") 
    if(cache)
    {   
        pwys <- .getResourceFromCache(pwys.name)
        if(!is.null(pwys)) return(pwys)
    }   
    pwys <- KEGGREST::keggList("pathway", org)
    message("download:")
    kgmls <- vapply(names(pwys), 
        function(pwy)
        {
            message(pwy)
            ind <- which(names(pwys) == pwy)
            if(ind == round(length(pwys) / 2)) Sys.sleep(10) 
            pwy <- sub("^path:", "", pwy)
            kgml <- KEGGREST::keggGet(pwy, "kgml")
            return(kgml)
        }, character(1))
    names(kgmls) <- sub("^path:", "", names(kgmls))

    if(is.null(out.dir))
    {
        pwys <- sapply(kgmls, parseKGML)
        .cacheResource(pwys, pwys.name)
        return(pwys)
    }

    message("write kgml files ...")
    if(!file.exists(out.dir)) dir.create(out.dir)
    for(n in names(kgmls)) 
        cat(kgmls[n], file=file.path(out.dir, paste(n, ".xml", sep="")))

    if(zip) 
    {
        message("zip ...")
        pattern <- file.path(out.dir, "*.xml")
        zip.file <- file.path(out.dir, paste(org, "zip", sep="."))
        system(paste("zip", zip.file, pattern))
        for(f in list.files(out.dir, pattern=pattern)) file.remove(f)
    }
}


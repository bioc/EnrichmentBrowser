############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-06-23 10:25:54
# 
# descr: configuration
# 
############################################################

.onLoad <- function(libname, pkgname) .setupConfig()

.setupConfig <- function()
{
    # (1) important colData, rowData, and result column names
    .ebrowser_config_cache[["PVAL.COL"]] <-
        .ebrowser_config_cache[["GSP.COL"]] <- "PVAL"
    .ebrowser_config_cache[["ADJP.COL"]] <- "ADJ.PVAL"
    .ebrowser_config_cache[["FC.COL"]] <- "FC"
    .ebrowser_config_cache[["PRB.COL"]] <- "PROBEID"
    .ebrowser_config_cache[["EZ.COL"]] <- "ENTREZID"
    .ebrowser_config_cache[["GN.COL"]] <- "GENENAME"
    .ebrowser_config_cache[["SYM.COL"]] <- "SYMBOL"
    .ebrowser_config_cache[["GRP.COL"]] <- "GROUP"
    .ebrowser_config_cache[["SMPL.COL"]] <- "SAMPLE"
    .ebrowser_config_cache[["BLK.COL"]] <- "BLOCK"
    .ebrowser_config_cache[["GS.COL"]] <- "GENE.SET"
    .ebrowser_config_cache[["PMID.COL"]] <- "PUBMED"

    # (2) URLs
    .ebrowser_config_cache[["NCBI.URL"]] <- "http://www.ncbi.nlm.nih.gov/"
    .ebrowser_config_cache[["PUBMED.URL"]] <-
        paste0(.ebrowser_config_cache[["NCBI.URL"]], "pubmed/")
    .ebrowser_config_cache[["GENE.URL"]] <-
        paste0(.ebrowser_config_cache[["NCBI.URL"]], "gene/")
    .ebrowser_config_cache[["KEGG.URL"]] <- "http://www.genome.jp/dbget-bin/"
    .ebrowser_config_cache[["KEGG.GENE.URL"]] <-
        paste0(.ebrowser_config_cache[["KEGG.URL"]], "www_bget?")
    .ebrowser_config_cache[["KEGG.SHOW.URL"]] <-
        paste0(.ebrowser_config_cache[["KEGG.URL"]], "show_pathway?")
    .ebrowser_config_cache[["GO.SHOW.URL"]] <- 
        "http://amigo.geneontology.org/amigo/term/"
    
    # (3) file paths
    .ebrowser_config_cache[["EBROWSER.HOME"]] <- 
        tools::R_user_dir("EnrichmentBrowser")
    .ebrowser_config_cache[["OUTDIR.DEFAULT"]] <- 
        file.path(.ebrowser_config_cache[["EBROWSER.HOME"]], "results")
    
    # (4) methodological defaults 
    .ebrowser_config_cache[["GS.MIN.SIZE"]] <- 5
    .ebrowser_config_cache[["GS.MAX.SIZE"]] <- 500

    # (5) output appearance
    .ebrowser_config_cache[["RESULT.TITLE"]] <- "Table of Results" 
    .ebrowser_config_cache[["NR.SHOW"]] <- 20
    .ebrowser_config_cache[["PLOT.WIDTH"]] <- 500
    .ebrowser_config_cache[["PLOT.HEIGHT"]] <- 500
        
    # (6) misc
    sbea.pkgs <- c("GSA", "PADOG", "globaltest", 
                    "GSVA", "EmpiricalBrownsMethod", "mgsa")
    names(sbea.pkgs) <- c("gsa", "padog", "globaltest", "gsva", "ebm", "mgsa")
    .ebrowser_config_cache[["SBEA.PKGS"]] <- sbea.pkgs
    nbea.pkgs <- c("PathNet", "DEGraph", "CePa", "GANPA", "NetGSA", "neaGUI")
    names(nbea.pkgs) <- c("pathnet", "degraph", "cepa", "ganpa", "netgsa", "nea")
    .ebrowser_config_cache[["NBEA.PKGS"]] <- nbea.pkgs
}


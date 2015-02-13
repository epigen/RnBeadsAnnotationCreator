########################################################################################################################
## createAnnotationPackage.hg38.R
## created: 2014-02-13
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for HG38.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

rnb.globals.hg38 <- function() {

	suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))

	## Supported chromosomes
	CHROMOSOMES <- c(1:22, "X", "Y")
	names(CHROMOSOMES) <- CHROMOSOMES
	assign('CHROMOSOMES', CHROMOSOMES, .globals)

	## Ensembl table for gene definitions
	assign('ENSEMBL.DATASET', "hsapiens_gene_ensembl", .globals)

	## Ensembl gene table attributes
	ENSEMBL.GENE.ATTRS <- c(
		"id" = "ensembl_gene_id",
		"chromosome" = "chromosome_name",
		"start" = "start_position",
		"end" = "end_position",
		"strand" = "strand",
		"symbol" = "mgi_symbol",
		"entrezID" = "entrezgene")
	assign('ENSEMBL.GENE.ATTRS' ENSEMBL.GENE.ATTRS, .globals)
}

########################################################################################################################

#' createAnnotationPackage.hg38
#' 
#' Helper function to create annotation package for genome assembly hg38
#' RnBeads annotation for that assembly
#' @param dest destination directory where the package should be generated
#' @return invisible \code{TRUE} if successful
#' @author Yassen Assenov, Fabian Mueller
createAnnotationPackage.hg38 <- function(dest=getwd()){
	invisible(TRUE)
}

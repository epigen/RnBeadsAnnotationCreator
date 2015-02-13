########################################################################################################################
## createAnnotationPackage.hg38.R
## created: 2014-02-13
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for HG38.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' createAnnotationPackage.hg38
#' 
#' Helper function to create annotation package for genome assembly hg38
#' RnBeads annotation for that assembly
#' @param dest destination directory where the package should be generated
#' @return invisible \code{TRUE} if successful
#' @author Fabian Mueller
#' @noRd
createAnnotationPackage.hg38 <- function(dest){

	suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))

	## Supported chromosomes
	CHROMOSOMES <- c(1:22, "X", "Y")
	names(CHROMOSOMES) <- CHROMOSOMES
	assign('CHROMOSOMES', CHROMOSOMES, .globals)

	## TODO: Download genes from Ensembl
	biomart.parameters <- list(
		database.name = "ensembl",
		dataset.name = "hsapiens_gene_ensembl",
		required.columns = c(
			"id" = "ensembl_gene_id",
			"chromosome" = "chromosome_name",
			"start" = "start_position",
			"end" = "end_position",
			"strand" = "strand",
			"symbol" = "mgi_symbol",
			"entrezID" = "entrezgene"))

	logger.start("Region Annotation")
	regions <- update.annot("regions", "region annotation", rnb.update.region.annotation, biomart.parameters)
	logger.completed()

}

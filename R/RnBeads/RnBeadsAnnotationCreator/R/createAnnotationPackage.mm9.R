########################################################################################################################
## createAnnotationPackage.mm9.R
## created: 2015-04-20
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for Mm9.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' createAnnotationPackage.mm9
#'
#' Helper function to create RnBeads annotation package for genome assembly mm9.
#'
#' @return None (invisible \code{NULL}).
#' @author Yassen Assenov
#' @noRd
createAnnotationPackage.mm9 <- function() {

	suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm9))

	## Genomic sequence and supported chromosomes
	GENOME <- "BSgenome.Mmusculus.UCSC.mm9"
	assign('GENOME', GENOME, .globals)
	CHROMOSOMES <- c(1:19, "X", "Y")
	names(CHROMOSOMES) <- paste0("chr", CHROMOSOMES)
	assign('CHROMOSOMES', CHROMOSOMES, .globals)
	rm(GENOME)

	## Define genomic regions
	biomart.parameters <- list(
		database.name = "ENSEMBL_MART_ENSEMBL",
		dataset.name = "mmusculus_gene_ensembl",
		required.columns = c(
			"id" = "ensembl_gene_id",
			"chromosome" = "chromosome_name",
			"start" = "start_position",
			"end" = "end_position",
			"strand" = "strand",
			"symbol" = "mgi_symbol",
			"entrezID" = "entrezgene"),
		host = "may2012.archive.ensembl.org")
	logger.start("Region Annotation")
	update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	rm(biomart.parameters)
	logger.completed()

	## Define genomic sites
	logger.start("Genomic Sites")
	update.annot("sites", "CpG annotation", rnb.update.sites)
	logger.completed()

	## Create all possible mappings from regions to sites
	logger.start("Mappings")
	update.annot("mappings", "mappings", rnb.create.mappings)
	logger.completed()

	## Export the annotation tables
	rnb.export.annotations.to.data.files()
}

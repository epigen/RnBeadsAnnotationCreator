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
	assignChromosomes(c(1:22, "X", "Y"))

	## Download SNP annotation
	logger.start("SNP Annotation")
	vcf.files <- paste0(DBSNP.FTP.BASE, "human_9606_b141_GRCh38/VCF/", "All.vcf.gz")
	snps <- update.annot("snps", "polymorphism information", rnb.update.dbsnp, ftp.files = vcf.files)
	logger.completed()

	## Download genes from Ensembl
	biomart.parameters <- list(
		database.name = "ENSEMBL_MART_ENSEMBL",
		dataset.name = "hsapiens_gene_ensembl",
		required.columns = c(
			"id" = "ensembl_gene_id",
			"chromosome" = "chromosome_name",
			"start" = "start_position",
			"end" = "end_position",
			"strand" = "strand",
			"symbol" = "hgnc_symbol",
			"entrezID" = "entrezgene"),
		host = "dec2014.archive.ensembl.org")

	logger.start("Region Annotation")
	regions <- update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	logger.completed()

	logger.start("Genomic Sites")
	sites <- update.annot("sites", "CpG annotation", rnb.update.sites, cpgislands = regions[["cpgislands"]])
	logger.completed()
}

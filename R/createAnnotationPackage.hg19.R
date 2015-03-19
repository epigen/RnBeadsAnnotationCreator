########################################################################################################################
## createAnnotationPackage.hg19.R
## created: 2014-02-13
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for HG19.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' createAnnotationPackage.hg19
#' 
#' Helper function to create annotation package for genome assembly hg19.
#'
#' @return invisible \code{TRUE} if successful
#' @author Fabian Mueller
#' @noRd
createAnnotationPackage.hg19 <- function() {

	suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))

	## Supported chromosomes
	assignChromosomes(c(1:22, "X", "Y"))

	## Download SNP annotation
	logger.start("SNP Annotation")
	vcf.files <- paste0(DBSNP.FTP.BASE, "human_9606_b141_GRCh37p13/VCF/", "00-All.vcf.gz")
	snps <- update.annot("snps", "polymorphism information", rnb.update.dbsnp, ftp.files = vcf.files)
	logger.info(paste("Using:", attr(snps, "version")))
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
		host = "feb2014.archive.ensembl.org")
	logger.start("Region Annotation")
	regions <- update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	logger.completed()

	## Define genomic sites
	logger.start("Genomic Sites")
	sites <- update.annot("sites", "CpG annotation", rnb.update.sites, cpgislands = regions[["cpgislands"]])
	logger.completed()

	## Create all possible mappings from regions to sites
	logger.start("Mappings")
	mappings <- update.annot("mappings", "mappings", rnb.create.mappings, regions = regions, sites = sites)
	logger.completed()

	

}

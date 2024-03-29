########################################################################################################################
## createAnnotationPackage.mm10.R
## created: 2014-04-20
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for Mm10.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' createAnnotationPackage.mm10
#'
#' Helper function to create RnBeads annotation package for genome assembly mm10.
#'
#' @return None (invisible \code{NULL}).
#' @author Yassen Assenov
#' @noRd
createAnnotationPackage.mm10 <- function() {

	suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm10))

	## Genomic sequence and supported chromosomes
	GENOME <- "BSgenome.Mmusculus.UCSC.mm10"
	assign('GENOME', GENOME, .globals)
	CHROMOSOMES <- c(1:19, "X", "Y")
	names(CHROMOSOMES) <- paste0("chr", CHROMOSOMES)
	assign('CHROMOSOMES', CHROMOSOMES, .globals)
	rm(GENOME)

	## Download SNP annotation
	logger.start("SNP Annotation")
	vcf.files <- gsub("^chr(.+)$", "vcf_chr_\\1.vcf.gz", names(CHROMOSOMES))
	vcf.files <- paste0(DBSNP.FTP.BASE.MM, "mouse_10090/VCF/", vcf.files)
	update.annot("snps", "polymorphism information", rnb.update.dbsnp, ftp.files = vcf.files)
	logger.info(paste("Using:", attr(.globals[['snps']], "version")))
	rm(vcf.files)
	logger.completed()

	## Define genomic regions
	biomart.parameters <- list(
		database.name = "ensembl",
		dataset.name = "mmusculus_gene_ensembl", db.version = 75,
		required.columns = c(
			"id" = "ensembl_gene_id",
			"chromosome" = "chromosome_name",
			"start" = "start_position",
			"end" = "end_position",
			"strand" = "strand",
			"symbol" = "mgi_symbol",
			"entrezID" = "entrezgene")) #"entrezgene_id" )) # for db.version=104
	logger.start("Region Annotation")
	update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	rm(biomart.parameters)
	logger.completed()

	## Define genomic sites
	logger.start("Genomic Sites")
	update.annot("sites", "CpG annotation", rnb.update.sites)
	logger.completed()
	
	## Define MethylationEPIC probe annotations
	logger.start("MouseMethylationBeadChip")
	table.columns <- rnb.get.illumina.annotation.columns("MOUSE")
	update.annot("probesMMBC", "MouseMethylationBeadChip annotation", rnb.update.probeMOUSE.annotation,
	             table.columns = table.columns)
	.globals[['sites']][["probesMMBC"]] <- .globals[['probesMMBC']][["probes"]]
	logger.completed()
	
	## Add annotation columns to the context probes, showing if they are covered by an assay
	logger.start("Updating Site Annotation with Probes")
	.globals[['sites']] <- rnb.update.site.annotation.with.probes(sites = .globals[['sites']],
                            query.probes = c("probesMMBC"),
                            platform.names = c("InfiniumMouseMBC"))
	logger.completed()

	## Create all possible mappings from regions to sites
	logger.start("Mappings")
	update.annot("mappings", "mappings", rnb.create.mappings)
	logger.completed()

	## Export the annotation tables
	rnb.export.annotations.to.data.files()
}

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
	suppressPackageStartupMessages(library(FDb.InfiniumMethylation.hg19))

	## Supported chromosomes
	assignChromosomes(c(1:22, "X", "Y"))

	## Download SNP annotation
	logger.start("SNP Annotation")
	vcf.files <- paste0(DBSNP.FTP.BASE, "human_9606_b141_GRCh37p13/VCF/", "00-All.vcf.gz")
	update.annot("snps", "polymorphism information", rnb.update.dbsnp, ftp.files = vcf.files)
	logger.info(paste("Using:", attr(.globals[['snps']], "version")))
	logger.completed()

	## Define genomic regions - tiling, genes (download from Ensembl), promoters, CpG islands (download from UCSC)
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
	update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	logger.completed()
	rm(biomart.parameters)

	## Define genomic sites
	logger.start("Genomic Sites")
	update.annot("sites", "CpG annotation", rnb.update.sites)
	logger.completed()

	## Define Infinium 27k probe annotations
	logger.start("Infinium 27k")
	ftp.table <- paste0(GEO.FTP.BASE, "GPL8490/GPL8490_HumanMethylation27_270596_v.1.2.csv.gz")
	table.columns <- c(
		"IlmnID" = "ID",
		"Name" = "Name",
		"IlmnStrand" = "IlmnStrand",
		"AddressA_ID" = "AddressA ID",
		"AlleleA_ProbeSeq" = "AlleleA Probe Sequence",
		"AddressB_ID" = "AddressB ID",
		"AlleleB_ProbeSeq" = "AlleleB Probe Sequence",
		"GenomeBuild" = "Genome Build",
		"Chr" = "Chromosome",
		"MapInfo" = "Location",
		"Ploidy" = "Ploidy",
		"Species" = "Species",
		"Source" = "Source",
		"SourceVersion" = "Source Version",
		"SourceStrand" = "Strand",
		"SourceSeq" = "Source Sequence",
		"TopGenomicSeq" = "Top Genomic Sequence",
		"Next_Base" = "Next Base",
		"Color_Channel" = "Color",
		"TSS_Coordinate" = "TSS Coordinate",
		"Gene_Strand" = "Gene Strand",
		"Gene_ID" = "Gene ID",
		"Symbol" = "Symbol",
		"Synonym" = "Synonym",
		"Accession" = "Accession",
		"GID" = "GID",
		"Annotation" = "Annotation",
		"Product" = "Product",
		"Distance_to_TSS" = "Distance to TSS",
		"CPG_ISLAND" = "CpG Island",
		"CPG_ISLAND_LOCATIONS" = "CpG Island Locations",
		"MIR_CPG_ISLAND" = "MIR CpG Island",
		"MIR_NAMES" = "MIR Names")
	update.annot("probes27", "Infinium 27K annotation", rnb.update.probe27k.annotation,
		ftp.table = ftp.table, table.columns = table.columns)
	.globals[['sites']][['probes27']] <- .globals[['probes27']][["probes"]]
	logger.completed()

	## Define Infinium 450k probe annotations
	logger.start("Infinium 450k")
	ftp.table <- paste0(GEO.FTP.BASE, "GPL13534/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
	table.columns <- c(
		"IlmnID" = "ID",
		"Name" = "Name",
		"AddressA_ID" = "AddressA_ID",
		"AlleleA_ProbeSeq" = "AlleleA Probe Sequence",
		"AddressB_ID" = "AddressB_ID",
		"AlleleB_ProbeSeq" = "AlleleB Probe Sequence",
		"Infinium_Design_Type" = "Design",
		"Next_Base" = "Next Base",
		"Color_Channel" = "Color",
		"Forward_Sequence" = "Forward Sequence",
		"Genome_Build" = "Genome Build",
		"CHR" = "Chromosome",
		"MAPINFO" = "Location",
		"SourceSeq" = "Source Sequence",
		"Chromosome_36" = "Chromosome.36",
		"Coordinate_36" = "Location.36",
		"Strand" = "Strand",
		"Probe_SNPs" = "Probe SNPs",
		"Probe_SNPs_10" = "Probe SNPs 10",
		"Random_Loci" = "Random",
		"Methyl27_Loci" = "HumanMethylation27",
		"UCSC_RefGene_Name" = "UCSC RefGene Name",
		"UCSC_RefGene_Accession" = "UCSC RefGene Accession",
		"UCSC_RefGene_Group" = "UCSC RefGene Group",
		"UCSC_CpG_Islands_Name" = "UCSC CpG Islands Name",
		"Relation_to_UCSC_CpG_Island" = "CGI Relation",
		"Phantom" = "Phantom",
		"DMR" = "DMR",
		"Enhancer" = "Enancer",
		"HMM_Island" = "HMM Island",
		"Regulatory_Feature_Name" = "Regulatory Feature Name",
		"Regulatory_Feature_Group" = "Regulatory Feature Group",
		"DHS" = "DHS")
	update.annot("probes450", "Infinium 450K annotation", rnb.update.probe450k.annotation,
		ftp.table = ftp.table, table.columns = table.columns)
	.globals[['sites']][['probes450']] <- .globals[['probes450']][["probes"]]
	logger.completed()

	## Create all possible mappings from regions to sites
	logger.start("Mappings")
	update.annot("mappings", "mappings", rnb.create.mappings)
	logger.completed()

	## Export the annotation tables
	rnb.export.annotations.to.data.files()
	return(invisible(TRUE))
}

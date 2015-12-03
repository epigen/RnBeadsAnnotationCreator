########################################################################################################################
## probes27.R
## created: 2015-04-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes the Infinium HumanMethylation27 probe definition tables by loading them from various sources.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.probe27k.annotation
#'
#' Creates probe annotation tables for Infinium 27k.
#'
#' @param ftp.table     FTP link to the probe annotation table in GEO.
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two items:
#'         \describe{
#'           \item{\code{"probes"}}{\code{GRangesList} instance containing probe annotation, one \code{GRanges} per
#'                chromosome.}
#'           \item{\code{"controls"}}{\code{data.frame} with control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probe27k.annotation <- function(ftp.table, table.columns) {

	## Download probe definition table from GEO
	probes27.geo <- rnb.load.probe.annotation.geo(ftp.table, table.columns, "HumanMethylation27")
	geo <- probes27.geo$probes

	## Validate the columns in the downloaded table
	if (anyDuplicated(geo[, "ID"]) != 0) {
		logger.error("Duplicated values for IllmnID")
	}
	rownames(geo) <- as.character(geo[, "ID"])
	geo[["Strand"]] <- rnb.fix.strand(geo[["Strand"]])
	geo[["IlmnStrand"]] <- rnb.fix.strand(geo[["IlmnStrand"]])
	if (!is.integer(geo[["AddressA"]])) {
		logger.error("Invalid values in column AddressA ID")
	}
	if (!is.integer(geo[["AddressB"]])) {
		logger.error("Invalid values in column AddressB ID")
	}
	geo[, "AlleleA Probe Sequence"] <- as.character(geo[["AlleleA Probe Sequence"]])
	geo[, "AlleleB Probe Sequence"] <- as.character(geo[["AlleleB Probe Sequence"]])
	if (!identical(levels(geo[["Color"]]), c("Grn", "Red"))) {
		logger.error("Invalid values in column Color")
	}
	levels(geo[["Color"]]) <- c("green", "red")

	## Obtain probe definition table from FDb.InfiniumMethylation.hg19
	probes27.bioc <- rnb.update.probe.annotation.methylumi("HumanMethylation27")
	methylumi <- probes27.bioc$probes # this table contains also 14 SNP probes
	methylumi <- methylumi[!is.na(methylumi[, "Chromosome"]), ]
	logger.status("Extracted probe definition table from FDb.InfiniumMethylation.hg19")

	probes.notmapped <- setdiff(rownames(geo), rownames(methylumi))
	if (length(probes.notmapped) != 0) {
		probes.notmapped <- paste(probes.notmapped, collapse = ", ")
		logger.warning(c("The following probes were not mapped to HG19:", probes.notmapped))
	}
	rm(probes.notmapped)

	## Construct data.frame combining the information from both sources
	probe.infos <- data.frame(
		"ID" = methylumi[, "ID"],
		"Chromosome" = methylumi[, "Chromosome"],
		"Location" = methylumi[, "Location"],
		"Strand" = geo[rownames(methylumi), "Strand"],
		"AddressA" = methylumi[, "AddressA"],
		"AddressB" = methylumi[, "AddressB"],
		"Design" = factor("I", levels = c("I", "II")),
		"IlmnStrand" = geo[rownames(methylumi), "IlmnStrand"],
		"AlleleA Probe Sequence" = geo[rownames(methylumi), "AlleleA Probe Sequence"],
		"AlleleB Probe Sequence" = geo[rownames(methylumi), "AlleleB Probe Sequence"],
		"Color" = geo[rownames(methylumi), "Color"],
		check.names = FALSE, stringsAsFactors = FALSE)
	rownames(probe.infos) <- probe.infos[, "ID"]
	probe.infos[is.na(probe.infos[["Strand"]]), "Strand"] <- "*"
	probe.infos[is.na(probe.infos[["IlmnStrand"]]), "IlmnStrand"] <- "*"
	i <- which(is.na(probe.infos[, "AddressB"]))
	if (length(i) != 0) {
		probe.infos[i, "AddressB"] <- probe.infos[i, "AddressA"]
	}
	logger.status("Combined both probe annotation tables")

	## Sort based on chromosome and position
	probe.infos <- probe.infos[with(probe.infos, order(Chromosome, Location)), ]

	## Add annotation for CpG density, GC content, sequence mismatches and SNPs
	probe.infos <- rnb.update.probe.annotation.cpg.context(probe.infos)
	probe.infos <- rnb.update.probe.annotation.msnps(probe.infos, .globals[['snps']])

	## Add data on cross-hybridization
	probe.infos[, "Cross-reactive"] <- rnb.update.probe.annotation.cr(probe.infos[, "ID"], "HumanMethylation27")

	## Convert to GRangesList
	probes.gr <- rnb.probe.infos.to.GRanges(probe.infos)

	## Validate and combine control probe annotations
	logger.start("Control Probes")
	geo <- probes27.geo$controls
	methylumi <- probes27.bioc$controls
	if (ncol(geo) != 4) {
		logger.error("Unexpected structure of control probes from GEO")
	}
	if (ncol(methylumi) != 4) {
		logger.error("Unexpected structure of control probes from Bioconductor")
	}
	colnames(geo) <- colnames(methylumi)
	geo$Type <- as.character(geo$Type)
	geo$Color_Channel <- as.character(geo$Color_Channel)
	if (anyDuplicated(geo$Address) != 0) {
		logger.error("Duplicated probe addresses in control probes from GEO")
	}
	if (anyDuplicated(methylumi$Address) != 0) {
		logger.error("Duplicated probe addresses in control probes from GEO")
	}
	geo <- geo[order(geo$Address), ]
	methylumi <- methylumi[order(methylumi$Address), ]
	if (!identical(geo$Address, methylumi$Address)) {
		logger.error("Incosistent probe addresses between GEO and Bioconductor")
	}
	if (identical(toupper(geo$Type), methylumi$Type) && identical(geo$Color_Channel, methylumi$Color_Channel)) {
		logger.info("GEO and Bioconductor agree on control probe annotation")
	} else {
		logger.warning("GEO and Bioconductor do NOT agree on control probe annotation; GEO will be taken")
	}
	logger.completed()

	return(list(probes = probes.gr, controls = geo))
}

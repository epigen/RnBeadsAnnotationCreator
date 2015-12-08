########################################################################################################################
## probesEPIC.R
## created: 2015-11-06
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes the Infinium MethylationEPIC probe definition tables by loading them from the Illumina website.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.probeEPIC.annotation
#'
#' Creates probe annotation tables for MethylationEPIC.
#'
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two items:
#'         \describe{
#'           \item{\code{"probes"}}{\code{GRangesList} instance containing probe annotations, one \code{GRanges} per
#'                chromosome.}
#'           \item{\code{"controls"}}{\code{data.frame} with control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probeEPIC.annotation <- function(table.columns) {

	ftp.table <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-manifest-file-csv.zip"

	## Download probe definition table from Illumina's web site
	destfile <- file.path(.globals[['DIR.PACKAGE']], "temp", "probeEPIC.zip")
	if (file.exists(destfile)) {
		logger.status(c("File", destfile, "already downloaded"))
	} else {
		if (download.file(ftp.table, destfile, quiet = TRUE, mode = "wb") != 0) {
			logger.error(c("Could not download", ftp.table))
		}
		logger.status(c("Downloaded", ftp.table))
	}

	## Unzip the downloaded archive
	dir.epic <- file.path(.globals[['DIR.PACKAGE']], "temp", "EPIC")
	if (file.exists(dir.epic)) {
		if (unlink(dir.epic, TRUE, TRUE) != 0) {
			logger.error(c("Could not remove existing temporary directory", dir.epic))
		}
	}
	if (!dir.create(dir.epic, FALSE)) {
		logger.error(c("Could not create temporary directory", dir.epic))
	}
	result <- unzip(destfile, exdir = dir.epic)
	if (length(result) != 1 || isTRUE(file.info(result)[, "isdir"])) {
		logger.error(c("Unexpected contents of", destfile))
	}

	## Identify methylation and control probe annotation tables
	txt <- scan(result, what = character(), sep = "\n", quiet = TRUE)
	assay.start <- grep("^\\[Assay\\]", txt)
	assay.start <- ifelse(length(assay.start) != 0, assay.start[1], 0L)
	controls.start <- grep("^\\[Controls\\]", txt)
	if (!(length(controls.start) == 1 && controls.start > assay.start + 1)) {
		logger.error("Missing or invalid [Controls] section")
	}
	rm(txt); invisible(gc())

	## Load the methylation and control probe annotation tables
	probe.infos <- read.csv(result, skip = assay.start, nrows = controls.start - assay.start - 2, check.names = FALSE)
	probe.infos <- probe.infos[, sapply(probe.infos, function(x) { !all(is.na(x)) })]
	if (!identical(colnames(probe.infos), names(table.columns))) {
		logger.error("Unexpected columns in the probe definition table")
	}
	colnames(probe.infos) <- table.columns[colnames(probe.infos)]
	control.probe.infos <- read.csv(result, header = FALSE, skip = controls.start, stringsAsFactors = FALSE)
	control.probe.infos <- control.probe.infos[, sapply(control.probe.infos, function(x) { !all(is.na(x)) })]
	logger.status("Loaded the probe definition tables from Illumina's web site")

	## Validate the control probes
	if (ncol(control.probe.infos) != 5) {
		logger.error("Unexpected number of columns in the control probe definition table")
	}
	colnames(control.probe.infos) <- INF.CONTROL.PROBE.TABLE.COLUMNS
	if (anyDuplicated(control.probe.infos$ID) != 0) {
		logger.error("Duplicated IDs in the control probe definition table")
	}
	if (!identical(sort(unique(control.probe.infos[, 2])), unname(RnBeads:::EPIC.CONTROL.TARGETS))) {
		logger.error("Unexpected values for Target in the control probe definition table")
	}
	control.probe.infos[, 2] <- factor(control.probe.infos[, 2], levels = unname(RnBeads:::EPIC.CONTROL.TARGETS))
	control.probe.infos[, 3] <- factor(RnBeads:::capitalize(tolower(control.probe.infos[, 3])))
	control.probe.infos[, 5] <- control.probe.infos[, 5] == "AVG"
	logger.status("Processed control probe annotation table")

	## Add information about CpG counts and GC content in the neighborhood, context, overlaps with SNPs
	probe.infos <- rnb.update.probe.annotation.cpg.context(probe.infos)
	probe.infos <- rnb.update.probe.annotation.msnps(probe.infos)

	## Add data on cross-hybridization
#	probe.infos[, "Cross-reactive"] <- rnb.update.probe.annotation.cr(probe.infos[, "ID"], "HumanMethylation450")

	## Convert to GRangesList
	probes.gr <- rnb.probe.infos.to.GRanges(probe.infos)

	return(list(probes = probes.gr, controls = control.probe.infos))
}

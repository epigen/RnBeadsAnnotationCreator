########################################################################################################################
## probes.R
## created: 2015-04-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Common functions for creating and updating Infinium probe annotation tables.
########################################################################################################################

rnb.load.probe.annotation.geo <- function(ftp.table, table.columns, platform = "HumanMethylation27k") {

	## Download probe definition table from GEO
	destfile <- file.path(.globals[['DATA.PACKAGE']], paste0(platform, ".csv.gz"))
	if (file.exists(destfile)) {
		logger.status(c("File", destfile, "already downloaded"))
	} else {
		if (download.file(ftp.table, destfile, quiet = TRUE, mode = "wb") != 0) {
			logger.error(c("Could not download", ftp.table))
		}
		logger.status(c("Downloaded", ftp.table))
	}
	txt <- scan(destfile, what = character(), sep = "\n", quiet = TRUE)
	assay.start <- grep("^\\[Assay\\]", txt)
	assay.start <- ifelse(length(assay.start) != 0, assay.start[1], 0L)
	controls.start <- grep("^\\[Controls\\]", txt)
	if (!(length(controls.start) == 1 && controls.start > assay.start + 1)) {
		logger.error("Missing or invalid [Controls] section")
	}
	# it is faster to re-read the file than to call read.csv from textConnection(txt)
	probe.infos <- read.csv(destfile, skip = assay.start, nrows = controls.start - assay.start - 2, check.names = FALSE)
	control.probe.infos <- read.csv(destfile, header = FALSE, skip = controls.start, check.names = FALSE)
	logger.status("Loaded the probe definition tables")

	## Validate probe.infos columns
	N <- length(table.columns)
	if (!identical(colnames(probe.infos)[1:N], names(table.columns))) {
		logger.error("Unexpected columns in the probe definition table")
	}

	probe.infos <- probe.infos[, sapply(probe.infos, function(x) { !all(is.na(x)) })]
	colnames(probe.infos) <- table.columns[colnames(probe.infos)]
	control.probe.infos <- control.probe.infos[, sapply(control.probe.infos, function(x) { !all(is.na(x)) })]
	return(list(probes = probe.infos, controls = control.probe.infos))
}

########################################################################################################################

#' rnb.update.probe.annotation.methylumi
#'
#' Creates probe annotation \code{data.frame}s using the data avaible in the package
#' \pkg{IlluminaHumanMethylation450k.db}.
#'
#' @param platform Assay name, must be one of \code{"HumanMethylation27k"} or \code{"HumanMethylation450k"}.
#' @return \code{list} of two elements:
#'         \describe{
#'           \item{\code{probes}}{Table of annotation for all methylation and SNP-based probes in the assay. It
#'                contains the following columns: \code{"Chromosome"}, \code{"Location"}, \code{"Random"},
#'                \code{"CGI Relation"} and \code{"CGI Location"}}
#'           \item{\code{controls}}{Table of annotation for all control probes in the assay, containing the
#'                following columns: \code{"Address"}, \code{"Type"}, \code{"Color_Channel"} and \code{"Name"}.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.methylumi <- function(platform = "HumanMethylation27k") {

	if (platform == "HumanMethylation27k") {
		probe.infos <- get27k() # mcols()[, "addressA", "addressB", "channel", "probeType", "platform", "percentGC"]
	} else { # platform == "HumanMethylation450k"
		probe.infos <- get450k() # mcols()[, "addressA", "addressB", "channel", "probeType", "platform", "percentGC"]
	}
	probe.infos <- data.frame(
		ID = names(probe.infos),
		Chromosome = as.character(seqnames(probe.infos)),
		Location = start(probe.infos),
		Context = as.factor(mcols(probe.infos)[, "probeType"]),
		Channel = as.factor(mcols(probe.infos)[, "channel"]),
		AddressA = as.integer(mcols(probe.infos)[, "addressA"]),
		AddressB = as.integer(mcols(probe.infos)[, "addressB"]),
		check.names = FALSE, stringsAsFactors = FALSE)
	rownames(probe.infos) <- probe.infos[, "ID"]

	## Set chromosome names
	chromosomes <- names(.globals[['CHROMOSOMES']])
	i <- which(!(probe.infos[, "Chromosome"] %in% chromosomes))
	if (length(i) != 0) {
		i <- paste(probe.infos[i, "ID"], collapse = ", ")
		logger.warning(c("The following probes are not located on unsupported chromosomes:", i))
	}
	probe.infos[, "Chromosome"] <- suppressWarnings(factor(probe.infos[, "Chromosome"], levels = chromosomes))

	## Sort based on chromosome and position
	probe.infos <- probe.infos[with(probe.infos, order(Chromosome, Location)), ]
	logger.status("Extracted probe annotation from the methylumi package")

	## Extract a table with control probe annotations
	if (platform == "HumanMethylation27k") {
		control.probe.infos <- hm27.controls # [, "Address", "Type", "Color_Channel", "Name"]
	} else { # platform == "HumanMethylation450k"
		control.probe.infos <- hm450.controls # [, "Address", "Type", "Color_Channel", "Name"]
	}
	rownames(control.probe.infos) <- control.probe.infos[, "Address"]
	logger.status("Extracted control probe annotation from the methylumi package")

	return(list("probes" = probe.infos, "controls" = control.probe.infos))
}

########################################################################################################################
## probes450.R
## created: 2015-04-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes the Infinium HumanMethylation450 probe definition tables by loading them from various sources.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.probe450k.annotation
#'
#' Creates probe annotation tables for Infinium 450k.
#'
#' @param ftp.table     FTP link to the probe annotation table in GEO.
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two items:
#'         \describe{
#'           \item{\code{"probes"}}{\code{GRangesList} instance containing probe annotations, one \code{GRanges} per
#'                chromosome.}
#'           \item{\code{"controls"}}{\code{data.frame} with control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probe450k.annotation <- function(ftp.table, table.columns) {

	logger.start("Loading Probe Definitions")
	geo <- rnb.update.probe.annotation.geo(ftp.table, table.columns)
	methylumi <- rnb.update.probe.annotation.methylumi("HumanMethylation450")
	compiled <- rnb.update.probe.annotation.compiled()
	logger.completed()

	controls450 <- rnb.update.controls450.combine(methylumi$controls, compiled, geo$controls)
	probes450 <- rnb.update.probes450.combine(methylumi$probes, geo$probes)
	logger.status("Updated information on probe sequence overlap with SNPs")
	return(list(probes = probes450, controls = controls450))
}

########################################################################################################################

#' rnb.update.probe.annotation.geo
#'
#' Creates probe annotation \code{data.frame}s using the Infinium platform definition in GEO.
#'
#' @param base.dir  Local directory to store the downloaded HumanMethylation450 probe definition table.
#' @return \code{list} of two elements:
#'         \describe{
#'           \item{\code{"probes"}}{Table of annotation for all methylation and SNP-based probes in the assay.}
#'           \item{\code{"controls"}}{Table of annotation for all control probes in the assay.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.geo <- function(ftp.table, table.columns) {

	## Download probe definition table from GEO
	result <- rnb.load.probe.annotation.geo(ftp.table, table.columns, "HumanMethylation450")
	probe.infos <- result$probes
	control.probe.infos <- result$controls
	rm(result)

	## Validate control.probe.infos columns and some of the values
	if (ncol(control.probe.infos) != 5 && is.integer(control.probe.infos[[1]])) {
		logger.error("Unexpected number of columns for control probes")
	}
	if (!setequal(RnBeads:::HM450.CONTROL.TARGETS, levels(control.probe.infos[[2]]))) {
		logger.error("Unexpected values for Target")
	}
	colnames(control.probe.infos) <- INF.CONTROL.PROBE.TABLE.COLUMNS
	rownames(control.probe.infos) <- control.probe.infos[, 1]
	control.probe.infos[, "Description"] <- as.character(control.probe.infos[, "Description"])
	control.probe.infos[, "AVG"] <- (control.probe.infos[, "AVG"] == "AVG")

	## Validate probe.infos columns and some of the values
	probe.infos <- rnb.probes.fix.infinium.columns(probe.infos)

	## Drop columns that are not used by RnBeads
	probe.infos <- probe.infos[, c("ID", "Design", "Color", "Random", "HumanMethylation27",
			"Genome Build", "Chromosome", "Location", "Strand", "CGI Relation",
			"AlleleA Probe Sequence", "AlleleB Probe Sequence", "AddressA", "AddressB")]
	logger.info("Dropped unused columns")

	return(list("probes" = probe.infos, "controls" = control.probe.infos))
}

########################################################################################################################

#' rnb.update.probe.annotation.compiled
#'
#' Creates a control probe annotation \code{data.frame} by loading and validating the control probe annotation table
#' compiled by Pavlo Lutsik.
#'
#' @return Control probe annotation table in the form of a \code{data.frame} with columns as specified in the compiled
#'         annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.compiled <- function() {
	filename <- system.file("extdata/control.meta.data.csv", package = "RnBeadsAnnotationCreator")
	if (!file.exists(filename)) {
		logger.error(c("Required file not found:", filename))
	}
	control.probe.infos <- read.csv(filename)
	logger.status(c("Loaded compiled table of control probe annotation table from", filename))

	## Validate column names in the CSV file
	expected.columns <- c("Address", "Purpose", "Description", "Evaluate.Green", "Evaluate.Red", "Expected.Intensity",
		"Control.Type", "TargetID", "Index")
	if (!identical(colnames(control.probe.infos), expected.columns)) {
		logger.error("Unexpected columns in the control probe annotation table")
	}
	## Adjust table values
	rownames(control.probe.infos) <- control.probe.infos[, "Address"]
	control.probe.infos[, "Description"] <- as.character(control.probe.infos[["Description"]])
	for (cname in c("Evaluate.Green", "Evaluate.Red")) {
		control.probe.infos[, cname] <- factor(as.character(control.probe.infos[, cname]), levels = c("-", "+"))
	}

	return(control.probe.infos)
}

########################################################################################################################

#' rnb.update.controls450.enrich
#'
#' Extends the given table of control probe annotations by adding (or replacing) four columns and filling them with
#' values depending on the values of some of the original columns.
#'
#' @param control.probe.infos Table of control probe annotation in the form of a \code{data.frame} containing at least
#'                            the following columns: \code{"Target"}, \code{"Color"} and \code{"Description"}.
#' @return Enriched probe annotation; the given \code{data.frame} with four added or replaced columns:
#'         \code{"Evaluate Green"}, \code{"Evaluate Red"}, \code{"Expected Intensity"} and \code{"Sample-dependent"}.
#' @author Yassen Assenov
#' @noRd
rnb.update.controls450.enrich <- function(control.probe.infos) {

	## Add columns Evaluate Green and Evaluate Red
	control.probe.infos[["Evaluate Green"]] <- as.character(NA)
	control.probe.infos[["Evaluate Red"]] <- as.character(NA)
	i <- grep("^DNP", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^Biotin", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.GREEN)
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.RED)
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^NEGATIVE", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	control.probe.infos[["Evaluate Green"]] <- factor(control.probe.infos[["Evaluate Green"]], levels = c("-", "+"))
	control.probe.infos[["Evaluate Red"]] <- factor(control.probe.infos[["Evaluate Red"]], levels = c("-", "+"))

	## Add column Expected Intensity
	control.probe.infos[["Expected Intensity"]] <- as.character(NA)
	i <- control.probe.infos[, "Target"] %in% c("NEGATIVE", "TARGET REMOVAL")
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- c("BISULFITE CONVERSION II", "SPECIFICITY II", "EXTENSION", "NON-POLYMORPHIC")
	i <- control.probe.infos[, "Target"] %in% i
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% paste("NORM", c("A", "C", "G", "T"), sep = "_")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((High)|(20K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\(Medium\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Medium"
	i <- grep("\\(Low\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Low"
	i <- grep("\\((Bkg)|(5K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^BS Conversion I C", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^BS Conversion I U", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^GT Mismatch.+\\(PM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^GT Mismatch.+\\(MM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	control.probe.infos[["Expected Intensity"]] <- factor(control.probe.infos[["Expected Intensity"]])

	## Add column Sample-dependent
	control.probe.infos[["Sample-dependent"]] <-
		!(control.probe.infos[["Target"]] %in% RnBeads:::CONTROL.TARGETS.SAMPLE.INDEPENDENT)

	return(control.probe.infos)
}

########################################################################################################################

#' rnb.update.controls450.combine
#'
#' Validates that the obtained control probe annotation tables are consistent and combines them.
#'
#' @param methylumi Table of control probe annotation as obtained from \pkg{methylumi}. The control probes are expected
#'                  to be sorted by identifier and contain the following columns: \code{"Address"}, \code{"Type"},
#'                  \code{"Color_Channel"} and \code{"Name"}.
#' @param compiled  Table of control probe annotation loaded from the table compiled by Lutsik. See
#'                  \code{\link{rnb.update.probe.annotation.compiled}} for more details.
#' @param geo       Table of control probe annotation downloaded from GEO using the \pkg{GEOquery} package.
#' @return The combined and enriched control probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.controls450.combine <- function(methylumi, compiled, geo) {

	## Validate methylumi match geo
	methylumi.match <- FALSE
	if (setequal(geo[["ID"]], methylumi[["Address"]])) {
		cpi <- methylumi[rownames(geo), ]
		methylumi.match <- identical(as.character(geo[["Target"]]), cpi[["Type"]]) &&
			identical(as.character(geo[["Color"]]), cpi[["Color_Channel"]]) &&
			identical(as.character(geo[["Description"]]), cpi[["Name"]])
		rm(cpi)
	}
	if (methylumi.match) {
		logger.info("Table of control probes matches the one in methylumi")
	} else {
		logger.warning("Table of control probes does NOT match the one in methylumi")
	}
	rm(methylumi.match)

	## Make descriptions consistent
	descr <- geo[, "Description"]
	descr <- sub("([^ ])\\(", "\\1 (", descr)
	descr <- sub("^BS Conversion I-", "BS Conversion I ", descr)
	geo[, "Description"] <- descr

	geo <- rnb.update.controls450.enrich(geo)

	## Validate geo matches the compiled table
	if (is.null(compiled)) {
		logger.warning("No compiled table of control probes provided; comparison is skipped")
	} else {
		if (setequal(rownames(geo), rownames(compiled))) {
			geo <- geo[rownames(compiled), ]
			if (identical(geo[, "Target"], compiled[, "TargetID"]) &&
				identical(geo[, "Evaluate Red"], compiled[, "Evaluate.Red"]) &&
				identical(geo[, "Evaluate Green"], compiled[, "Evaluate.Green"]) &&
				identical(geo[, "Sample-dependent"], compiled[, "Control.Type"] == "Sample-Dependent") &&
				identical(geo[, "Expected Intensity"], compiled[, "Expected.Intensity"])) {
				logger.info("Table of control probes matches the one in the compiled table")
				geo[, "Index"] <- compiled[, "Index"]
				logger.info("Copied column Index to the table of control probes")
			} else {
				logger.warning("Table of control probes does NOT match the one in the compiled table")
			}
		} else {
			logger.warning("Control probe identifiers do NOT match the ones in the compiled table")
		}
	}
	return(geo)
}

########################################################################################################################

#' rnb.update.probes450.combine
#'
#' Validates that the obtained probe annotation tables are consistent and combines them.
#'
#' @param methylumi Probe annotation table obtained from the \pkg{methylumi} package, in the form of a
#'                  \code{data.frame}.
#' @param geo       Probe annotation table downloaded from GEO using the \pkg{GEOquery} package, in the form of a
#'                  \code{data.frame}.
#' @return The combined and enriched probe annotation table.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.probes450.combine <- function(methylumi, geo) {

	## Copy chromosome location of rs probes from methylumi to geo
	probes.rs <- mapply(grep, x = list(methylumi = rownames(methylumi), geo = rownames(geo)), pattern = "^rs",
		value = TRUE, SIMPLIFY = FALSE)
	if (setequal(probes.rs$methylumi, probes.rs$geo)) {
		probes.rs <- probes.rs$geo
		if (length(probes.rs) != 0) {
			location.rs <- cbind(geo = geo[probes.rs, "Location"], methylumi = methylumi[probes.rs, "Location"])
			i <- which(is.na(location.rs[, "geo"]) & (!is.na(location.rs[, "methylumi"])))
			if (length(i) != 0) {
				geo[probes.rs[i], "Chromosome"] <- methylumi[probes.rs[i], "Chromosome"]
				geo[probes.rs[i], "Location"] <- methylumi[probes.rs[i], "Location"]
				geo <- geo[with(geo, order(Chromosome, Location)), ]
				logger.status(c("Copied information for ", length(i), "rs probes from methylumu to geo"))
			}
			rm(location.rs, i)
		}
	} else {
		logger.warning("Mismatch between rs probe annotations from methylumi and geo; methylumi annotation is ignored")
	}

	## Validate geo matches methylumi
	if (setequal(geo[["ID"]], methylumi[["ID"]])) {
		probe.infos.methylumi <- methylumi[geo[["ID"]], ]
		col.notcovered <- which(!(colnames(probe.infos.methylumi) %in% colnames(geo)))
		if (length(col.notcovered) != 0) {
			msg <- paste(colnames(probe.infos.methylumi)[col.notcovered], collapse = ", ")
			logger.warning(c("The following column(s) from methylumi's probe annotation are ignored:", msg))
			probe.infos.methylumi <- probe.infos.methylumi[, -col.notcovered]
			rm(msg)
		}
		if (identical(geo[, colnames(probe.infos.methylumi)], probe.infos.methylumi)) {
			logger.info("Probe annotations match the data in the table from methylumi")
		} else {
			logger.warning("Probe annotations do NOT match the data in the table from methylumi")
		}
		rm(probe.infos.methylumi, col.notcovered)
	} else {
		logger.warning(c("Probe identifiers do NOT match the ones in the compiled table"))
	}

	geo[["Mismatches A"]] <- as.integer(NA)
	geo[["Mismatches B"]] <- as.integer(NA)
	if (!("Context" %in% colnames(geo))) {
		geo[["Context"]] <- as.character(NA)
	}
	geo <- geo[, c("ID", "Design", "Color", "Context", "Random", "HumanMethylation27",
		"Mismatches A", "Mismatches B", "Genome Build", "Chromosome", "Location", "Strand", "CGI Relation",
		"AlleleA Probe Sequence", "AlleleB Probe Sequence", "AddressA", "AddressB")]
	i <- which(is.na(geo[, "AddressB"]))
	if (length(i) != 0) {
		geo[i, "AddressB"] <- geo[i, "AddressA"]
	}

	## Add information about CpG counts and GC content in the neighborhood, context, overlaps with SNPs
	geo <- rnb.update.probe.annotation.cpg.context(geo)
	geo <- rnb.update.probe.annotation.msnps(geo)

	## Add data on cross-hybridization
	geo[, "Cross-reactive"] <- rnb.update.probe.annotation.cr(geo[, "ID"], "HumanMethylation450")

	## Convert to GRangesList
	probes.gr <- rnb.probe.infos.to.GRanges(geo)

	return(probes.gr)
}

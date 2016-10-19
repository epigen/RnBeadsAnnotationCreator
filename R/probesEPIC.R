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

	ftp.table <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v1-0-b2-manifest-file-csv.zip"

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
	probe.infos <- read.csv(result, skip = assay.start, nrows = controls.start - assay.start - 2, check.names = FALSE,
		stringsAsFactors = FALSE)
	probe.infos <- probe.infos[, sapply(probe.infos, function(x) { !all(is.na(x)) })]
	if (!identical(colnames(probe.infos), names(table.columns))) {
		logger.error("Unexpected columns in the probe definition table")
	}
	colnames(probe.infos) <- table.columns[colnames(probe.infos)]
	control.probe.infos <- read.csv(result, header = FALSE, skip = controls.start, stringsAsFactors = FALSE)
	control.probe.infos <- control.probe.infos[, sapply(control.probe.infos, function(x) { !all(is.na(x)) })]
	logger.status("Loaded the probe definition tables from Illumina's web site")

	## Validate probe.infos columns and some of the values
	probe.infos[, "Design"] <- as.factor(probe.infos[, "Design"])
	probe.infos[, "Color"] <- as.factor(probe.infos[, "Color"])
	probe.infos[, "Next Base"] <- as.factor(probe.infos[, "Next Base"])
	probe.infos[, "CGI Relation"] <- as.factor(probe.infos[, "CGI Relation"])
	probe.infos[, "DMR"] <- as.factor(probe.infos[, "DMR"])
	probe.infos[, "Enhancer"] <- as.factor(probe.infos[, "Enhancer"])
	probe.infos[, "Regulatory Feature Group"] <- as.factor(probe.infos[, "Regulatory Feature Group"])
	probe.infos <- rnb.probes.fix.infinium.columns(probe.infos)

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
	control.probe.infos <- rnb.update.controlsEPIC.enrich(control.probe.infos)
	logger.status("Processed control probe annotation table")

	## Add information about CpG counts and GC content in the neighborhood, context, overlaps with SNPs
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-1.RDS"))
	probe.infos <- rnb.update.probe.annotation.guess.strand(probe.infos)
	i <- which(probe.infos[, "Strand"] != probe.infos[, "Guessed Strand"])
	if (length(i) != 0) {
		i <- table(substr(rownames(probe.infos)[i], 1, 2))
		i <- paste0(names(i), " (", i, " probes)")
		logger.warning(paste("The Strand info for the following probe types might be wrong:", i))
	}
	rm(i)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-2.RDS"))
	probe.infos <- rnb.update.probesEPIC.snps(probe.infos)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-3.RDS"))
	probe.infos <- rnb.update.probe.annotation.cpg.context(probe.infos)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-4.RDS"))
	probe.infos <- rnb.update.probe.annotation.snps(probe.infos)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-5.RDS"))

	## Add data on cross-hybridization
	probe.infos[, "Cross-reactive"] <- rnb.update.probe.annotation.cr(probe.infos[, "ID"], "MethylationEPIC")

	## Convert to GRangesList
	probes.gr <- rnb.probe.infos.to.GRanges(probe.infos)

	return(list(probes = probes.gr, controls = control.probe.infos))
}

########################################################################################################################

#' rnb.update.controlsEPIC.enrich
#'
#' Extends the given table of control probe annotations by adding (or replacing) four columns and filling them with
#' values depending on the values of some of the original columns.
#'
#' @param control.probe.infos Table of control probe annotation in the form of a \code{data.frame} containing at least
#'                            the following columns: \code{"Target"}, \code{"Color"} and \code{"Description"}.
#' @return Enriched probe annotation; the given \code{data.frame} with four added or replaced columns:
#'         \code{"Evaluate Green"}, \code{"Evaluate Red"}, \code{"Expected Intensity"} and \code{"Sample-dependent"}.
#' @author Pavlo Lutsik
#' @noRd
rnb.update.controlsEPIC.enrich <- function(control.probe.infos) {

	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.COLORS.GREEN <- c("Black", "Blue", "Cyan", "Green", "Lime", "Limegreen", "Skyblue")

	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.COLORS.RED <- c("Gold", "Orange", "Purple", "Red", "Tomato", "Yellow")

	## Add columns Evaluate Green and Evaluate Red
	control.probe.infos[["Evaluate Green"]] <- "-"
	control.probe.infos[["Evaluate Red"]] <- "-"
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
	i <- control.probe.infos[, "Target"] %in% c("NEGATIVE", "TARGET REMOVAL", "RESTORATION")
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- c("BISULFITE CONVERSION II", "SPECIFICITY II", "EXTENSION", "NON-POLYMORPHIC")
	i <- control.probe.infos[, "Target"] %in% i
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% paste("NORM", c("A", "C", "G", "T"), sep = "_")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((High)|(20K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((Medium)|(5K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Medium"
	i <- grep("\\(Low\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Low"

	#i <- grep("\\((Bkg)|(5K)\\)$", control.probe.infos[, "Description"])
	i <- grep("\\(Bkg\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^BS Conversion I[- ]C", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^BS Conversion I[- ]U", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^GT Mismatch.+\\(PM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^GT Mismatch.+\\(MM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	control.probe.infos[["Expected Intensity"]] <- factor(control.probe.infos[["Expected Intensity"]])

	## Add column Sample-dependent
	control.probe.infos[["Sample-dependent"]] <-
		!(control.probe.infos[["Target"]] %in% RnBeads:::CONTROL.TARGETS.SAMPLE.INDEPENDENT)

	## Add column Index
	control.probe.infos[["Index"]][order(control.probe.infos$Target)] <-
		unlist(sapply(sort(unique(control.probe.infos$Target)), function(target) {
				1:length(which(control.probe.infos$Target==target))
			}
		))

	return(control.probe.infos)
}

########################################################################################################################

#' rnb.update.probe.annotation.guess.strand
#' 
#' Updates the MethylationEPIC probe annotation table by adding a column named \code{Guessed Strand} and setting it to
#' the best guess based on the probe sequence and design type.
#' 
#' @param probe.infos Probe annotation table for MethylationEPIC in the form of a \code{data.frame}.
#' @return The updated probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.guess.strand <- function(probe.infos) {
	genome.data <- rnb.genome.data()
	guessed <- data.frame(
		"Guessed Strand" = factor(rep("*", nrow(probe.infos)), levels = c("+", "-", "*")),
		"Mismatches A" = 0L, "Mismatches B" = 0L, check.names = FALSE)

	for (chromosome in names(.globals[['CHROMOSOMES']])) {
		chrom.sequence <- genome.data[[chromosome]]
		for (pr.design in c("I", "II")) {
			i <- which(probe.infos[["Chromosome"]] == chromosome & probe.infos[["Design"]] == pr.design)
			if (length(i) != 0) {
				alleles.A <- as.character(probe.infos[i, "AlleleA Probe Sequence"])
				if (pr.design == "I") {
					alleles.B <- as.character(probe.infos[i, "AlleleB Probe Sequence"])
				} else {
					alleles.B <- NULL
				}
				loci <- probe.infos[i, "Location"]
				guessed[i, ] <- rnb.seq.guess.strands(chrom.sequence, loci, pr.design, alleles.A, alleles.B)
				rm(alleles.A, alleles.B, loci)
			}
		}
	}
	rm(chromosome, chrom.sequence, pr.design, i)

	i <- which(probe.infos$Context == "Other")
	guessed[i, "Guessed Strand"] <- "*"
	guessed[i, "Mismatches A"] <- NA
	guessed[i, "Mismatches B"] <- NA
	cbind(probe.infos, guessed)
}

########################################################################################################################

#' rnb.update.probesEPIC.snps
#' 
#' Sets chromosome and location information for the SNP probes in MethylationEPIC array, copying it from Infinium 450k.
#' 
#' @param probe.infos Probe annotation table for MethylationEPIC in the form of a \code{data.frame}.
#' @return The updated probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probesEPIC.snps <- function(probe.infos) {
	## Identify the SNP probes in Infinium 450k
	probes450.rs <- rnb.annotation2data.frame(.globals$sites$probes450)
	probes450.rs <- probes450.rs[grep("^rs", rownames(probes450.rs)), ]
	extended.ids.450k <- paste(probes450.rs$ID, probes450.rs$Design, probes450.rs$Color, sep = ".")
	
	## Map the SNP probes in EPIC to the corresponding SNP probes in Infinium 450k
	i <- grep("^rs", probe.infos$ID)
	extended.ids.epic <- paste(probe.infos$ID[i], probe.infos$Design[i], probe.infos$Color[i], sep = ".")
	if (!all(extended.ids.epic %in% extended.ids.450k)) {
		logger.error("Not all SNP probes in EPIC found in 450k")
	}
	i450k <- unname(sapply(extended.ids.epic, function(x) { which(extended.ids.450k == x) }))
	
	## Transfer the information about the SNP probes to EPIC
	probe.infos[i, "Genome Build"] <- 37L
	probe.infos[i, "Chromosome"] <- probes450.rs[i450k, "Chromosome"]
	probe.infos[i, "Location"] <- probes450.rs[i450k, "Start"]
	probe.infos[i, "Strand"] <- probes450.rs[i450k, "Strand"]
	probe.infos[i, "HumanMethylation450"] <- TRUE
	probe.infos
}

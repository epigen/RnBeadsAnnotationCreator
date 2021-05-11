########################################################################################################################
## probesMOUSE.R
## created: 2015-11-06
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes the Infinium MethylationMOUSE probe definition tables by loading them from the Illumina website.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.probeMOUSE.annotation
#'
#' Creates probe annotation tables for MethylationMOUSE.
#'
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two items:
#'         \describe{
#'           \item{\code{"probes"}}{\code{GRangesList} instance containing probe annotations, one \code{GRanges} per
#'                chromosome.}
#'           \item{\code{"controls"}}{\code{data.frame} with control probe annotation.}
#'           \item{\code{"flagged"}}{\code{GRangesList} with MFG change flagged probes.}
#'         }
#' @author Maximilian Schoenung
#' @noRd
rnb.update.probeMOUSE.annotation <- function(table.columns) {

  man.file <- "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/Infinium%20Mouse%20Methylation%20v1.0%20A1%20GS%20Manifest%20File.csv"

  ## create temporary directory
  dir.mouse <- file.path(.globals[['DIR.PACKAGE']], "temp", "MOUSE")
  if (file.exists(dir.mouse)) {
    if (unlink(dir.mouse, TRUE, TRUE) != 0) {
      logger.error(c("Could not remove existing temporary directory", dir.mouse))
    }
  }
  if (!dir.create(dir.mouse, FALSE,recursive = T)) {
    logger.error(c("Could not create temporary directory", dir.mouse))
  }
  
	## Download probe definition table from Illumina's web site
	destfile <- file.path(dir.mouse, "probeMOUSE.csv")
	if (file.exists(destfile)) {
		logger.status(c("File", destfile, "already downloaded"))
	} else {
		if (download.file(man.file, destfile, quiet = TRUE, mode = "wb") != 0) {
			logger.error(c("Could not download", man.file))
		}
		logger.status(c("Downloaded", man.file))
	}

	## Assign destfile to result
	result <- destfile

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
	#probe.infos2=probe.infos
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
	levels(probe.infos[, "Color"]) <- c("Both", "Grn", "Red")
	probe.infos[, "Next Base"] <- as.factor(probe.infos[, "Next Base"])
	probe.infos[, "Chromosome"] <- paste0("chr",probe.infos[, "Chromosome"])
	probe.infos[, "Chromosome"] <- factor(as.character(probe.infos[, "Chromosome"]),
	                                      levels=paste0("chr",c(0,1:19,"MT","X","Y")))
	probe.infos[, "Strand"] <- factor(probe.infos[, "Strand"],
	                                      labels=c("*","+","-"))
	probe.infos[,"Context"] <- factor(toupper(substr(probe.infos[,"Name"], start = 1, stop = 2)))
  
	
	## Validate the control probes
	if (ncol(control.probe.infos) != 4) {
		logger.error("Unexpected number of columns in the control probe definition table")
	}
	colnames(control.probe.infos) <- Mouse.CONTROL.PROBE.TABLE.COLUMNS
	if (anyDuplicated(control.probe.infos$ID) != 0) {
		logger.error("Duplicated IDs in the control probe definition table")
	}
	if (!identical(sort(unique(control.probe.infos[, 2])), unname(RnBeads:::EPIC.CONTROL.TARGETS))) {
		logger.error("Unexpected values for Target in the control probe definition table")
	}
	control.probe.infos[, 2] <- factor(control.probe.infos[, 2], levels = unname(RnBeads:::EPIC.CONTROL.TARGETS))
	control.probe.infos[, 3] <- factor(RnBeads:::capitalize(tolower(control.probe.infos[, 3])))
  #control.probe.infos[, 5] <- control.probe.infos[, 5] == "AVG"
	control.probe.infos <- rnb.update.controlsMOUSE.enrich(control.probe.infos)
	logger.status("Processed control probe annotation table")

	## Add information about CpG counts and GC content in the neighborhood, context, overlaps with SNPs
	#probe.infos[, "CGI Relation"] <- NA
	for(con in c("N Shelf","N Shore","CpG Island","S Shelf","S Shore")){
	  probe.infos[,con][probe.infos[, con]==""] <- NA
	}
	probes.shaped <- probe.infos[,c("Name","Chromosome","Location","Strand","Design","CpG Island","N Shelf","N Shore","S Shelf","S Shore",
	                                "AddressA","AddressB","Color","MFG Change Flagged","ID","Context")]
	probes.flagged <- probes.shaped[probes.shaped[,"MFG Change Flagged"]==TRUE,]
    #probes.flagged$Chromosome<-factor(as.character(probes.flagged$Chromosome), levels=setdiff(levels(probes.flagged$Chromosome), "chr0"))
	probes.shaped <- probes.shaped[probes.shaped[,"MFG Change Flagged"]==FALSE,]
    probes.shaped$Chromosome<-factor(as.character(probes.shaped$Chromosome), levels=setdiff(levels(probes.shaped$Chromosome), "chr0"))
    
    probes.shaped <- rnb.update.probe.annotation.snps(probes.shaped)
    
	# sort and GRanges List
	probes.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = probes.shaped,
	                                                         keep.extra.columns = T,
	                                                         seqnames.field = "Chromosome",
	                                                         start.field = "Location",
	                                                         end.field = "Location",
	                                                         strand.field = "Strand"
	                                                          )
	#probes.ranges <- IRanges::resize(probes.ranges,2)
	width(probes.ranges) <- width(probes.ranges)+1
	gr <- sort(probes.ranges, ignore.strand=TRUE)
	names(gr) <- gr$ID
	probes.gr <- GenomicRanges::split(gr, seqnames(gr))
	probes.gr <- probes.gr[names(.globals[['CHROMOSOMES']])]

	# sort and GRanges List for flagged
	flagged.ranges <- GenomicRanges::makeGRangesFromDataFrame(df = probes.flagged,
	                                                         keep.extra.columns = T,
	                                                         seqnames.field = "Chromosome",
	                                                         start.field = "Location",
	                                                         end.field = "Location",
	                                                         strand.field = "Strand"
	)
	width(flagged.ranges) <- width(flagged.ranges)+1
	fl <- sort(flagged.ranges)
	names(fl) <- fl$ID
	probes.fl <- GenomicRanges::split(fl, seqnames(fl))
	
	return(list(probes = probes.gr, controls = control.probe.infos, flagged=probes.fl))
}

########################################################################################################################

#' rnb.update.controlsMOUSE.enrich
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
rnb.update.controlsMOUSE.enrich <- function(control.probe.infos) {

	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.DES.GREEN <- c("ctl-BISULFITE-CONVERSION-I-140U_MUS","ctl-BISULFITE-CONVERSION-I-140M_MUS",
	                       "ctl-BISULFITE-CONVERSION-I-303U_MUS","ctl-BISULFITE-CONVERSION-I-303M_MUS",
	                       "ctl-BISULFITE-CONVERSION-I-2U_HSA","ctl-BISULFITE-CONVERSION-I-2M_HSA",
	                       "ctl-BISULFITE-CONVERSION-I-3U_HSA","ctl-BISULFITE-CONVERSION-I-3M_HSA",
	                       "ctl-EXTENSION-Extension-C","ctl-EXTENSION-Extension-G",
	                       "ctl-NON-POLYMORPHIC-NP-C_MUS","ctl-NON-POLYMORPHIC-NP-C_HSA",
	                       "ctl-NON-POLYMORPHIC-NP-G_MUS","ctl-NON-POLYMORPHIC-NP-G_HSA",
	                       "ctl-RESTORATION-Restore","ctl-SPECIFICITY-I-37","ctl-SPECIFICITY-I-39",
	                       "ctl-SPECIFICITY-I-GT-Mismatch-2-MM","ctl-SPECIFICITY-I-GT-Mismatch-2-PM"
	                       )

	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.DES.RED <- c("ctl-BISULFITE-CONVERSION-I-317U_MUS","ctl-BISULFITE-CONVERSION-I-317M_MUS",
	                     "ctl-BISULFITE-CONVERSION-I-318U_MUS","ctl-BISULFITE-CONVERSION-I-318M_MUS",
	                     "ctl-BISULFITE-CONVERSION-I-330U_MUS","ctl-BISULFITE-CONVERSION-I-330M_MUS",
	                     "ctl-BISULFITE-CONVERSION-I-4U_HSA","ctl-BISULFITE-CONVERSION-I-4M_HSA",
	                     "ctl-BISULFITE-CONVERSION-I-5U_HSA","ctl-BISULFITE-CONVERSION-I-5M_HSA",
	                     "ctl-BISULFITE-CONVERSION-I-6U_HSA","ctl-BISULFITE-CONVERSION-I-6M_HSA",
	                     "ctl-EXTENSION-Extension-A","ctl-EXTENSION-Extension-T",
	                     "ctl-NON-POLYMORPHIC-NP-T_MUS","ctl-NON-POLYMORPHIC-NP-T_HSA",
	                     "ctl-NON-POLYMORPHIC-NP-A_MUS","ctl-NON-POLYMORPHIC-NP-A_HSA",
	                     "ctl-SPECIFICITY-I-34",
	                     "ctl-SPECIFICITY-I-GT-Mismatch-4-MM","ctl-SPECIFICITY-I-GT-Mismatch-4-PM",
	                     "ctl-SPECIFICITY-I-GT-Mismatch-5-MM","ctl-SPECIFICITY-I-GT-Mismatch-5-PM"
	)

	## Add columns Evaluate Green and Evaluate Red
	control.probe.infos[["Evaluate Green"]] <- "-"
	control.probe.infos[["Evaluate Red"]] <- "-"
	i <- grep("DNP", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("Biotin", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Description"]] %in% CONTROL.DES.GREEN)
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Description"]] %in% CONTROL.DES.RED)
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("BISULFITE-CONVERSION-II",control.probe.infos[["Description"]])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("NEGATIVE", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("HYBRIDIZATION", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- grep("NORM_A", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("NORM_T", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("NORM_C", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- grep("NORM_G", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- grep("SPECIFICITY II", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("TARGET REMOVAL", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	control.probe.infos[["Evaluate Green"]] <- factor(control.probe.infos[["Evaluate Green"]], levels = c("-", "+"))
	control.probe.infos[["Evaluate Red"]] <- factor(control.probe.infos[["Evaluate Red"]], levels = c("-", "+"))

	## Add column Expected Intensity
	control.probe.infos[["Expected Intensity"]] <- as.character(NA)
	i <- control.probe.infos[, "Target"] %in% c("NEGATIVE", "RESTORATION")
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- control.probe.infos[, "Target"] %in% c("TARGET REMOVAL")
	control.probe.infos[i, "Expected Intensity"] <- "Low"
	i <- grep("U_", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("M_", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- c("BISULFITE CONVERSION II", "SPECIFICITY II", "EXTENSION", "NON-POLYMORPHIC")
	i <- control.probe.infos[, "Target"] %in% i
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% paste("NORM", c("A", "C", "G", "T"), sep = "_")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% c("SPECIFICITY I")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("-MM", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("Bkg", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("High", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("Medium", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Medium"
	i <- grep("Low", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Low"
	
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


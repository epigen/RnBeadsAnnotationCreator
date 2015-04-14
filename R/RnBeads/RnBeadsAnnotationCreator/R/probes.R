########################################################################################################################
## probes.R
## created: 2015-04-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Common functions for creating and updating Infinium probe annotation tables.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.load.probe.annotation.geo
#'
#' Loads the specified probe annotation table from GEO.
#'
#' @param ftp.table     FTP link to the probe annotation table in GEO.
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two \code{data.frame}s:
#'         \describe{
#'           \item{\code{"probes"}}{Probe annotation; column names are the ones specified in \code{table.columns}.}
#'           \item{\code{"controls"}}{Control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.load.probe.annotation.geo <- function(ftp.table, table.columns, platform = "HumanMethylation27k") {

	## Download probe definition table from GEO
	destfile <- file.path(.globals[['DIR.PACKAGE']], "temp", paste0(platform, ".csv.gz"))
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
	logger.status("Loaded the probe definition tables from GEO")

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

	pltrf <- ifelse(platform == "HumanMethylation27k", "HM27", "HM450")
	probe.infos <- suppressMessages(getPlatform(platform = pltrf, genome = .globals[['assembly']]))
	# mcols()[, "addressA", "addressB", "channel", "probeType", "platform", "percentGC"]
	probe.infos <- data.frame(
		ID = names(probe.infos),
		Chromosome = as.character(seqnames(probe.infos)),
		Location = start(probe.infos),
		Context = mcols(probe.infos)[, "probeType"],
		Channel = mcols(probe.infos)[, "channel"],
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

########################################################################################################################

#' rnb.update.probe.annotation.cpg.context
#'
#' Enriches the given probe annotation with CpG density, GC content and probe context.
#'
#' @param probe.infos Probe annotation table in the form of a \code{data.frame} containing at least the following
#'                    columns: \code{"ID"}, \code{"Chromosome"} and \code{"Location"}.
#' @return The modified probe annotation table. The resulting table contains the following columns added to or updated
#'         in the original one: \code{"CpG"}, \code{"GC"} and \code{"Context"}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.cpg.context <- function(probe.infos) {

	locations <- tapply(probe.infos[, "Location"], probe.infos[, "Chromosome"], identity, simplify = FALSE)

	genome.data <- get.genome.data()
	calculate.cg <- function(chrom, loci, l.neighborhood) {
		starts <- loci - l.neighborhood / 2L + 1L
		ends <- loci + l.neighborhood / 2L - 1L
		chrom.sequence <- genome.data[[chrom]]
		chrom.regions <- suppressWarnings(Views(chrom.sequence, start = starts, end = ends))

		result <- data.frame(
			"CpG" = dinucleotideFrequency(chrom.regions)[, "CG"],
			"GC" = as.integer(rowSums(letterFrequency(chrom.regions, c("C", "G")))))

		## Define probe context based on the targeted dinucleotide - CG, CC, CAG, CAH, CTG, CTH, Other
		p.contexts <- factor(rep("CG", length(loci)), levels = c("CG", "CC", "CAG", "CAH", "CTG", "CTH", "Other"))
		bases1 <- suppressWarnings(as.character(Views(chrom.sequence, loci, loci)))
		loci <- loci + 1L
		bases2 <- suppressWarnings(as.character(Views(chrom.sequence, loci, loci)))
		p.contexts[which(bases1 != "C")] <- "Other"
		i.noncg <- which(bases1 == "C" & bases2 != "G")
		if (length(i.noncg) != 0) {
			bases2 <- bases2[i.noncg]
			loci <- loci[i.noncg] + 1L
			bases3 <- suppressWarnings(as.character(Views(chrom.sequence, loci, loci)))
			p.contexts[i.noncg[bases2 == "C"]] <- "CC"
			p.contexts[i.noncg[bases2 == "A" & bases3 == "G"]] <- "CAG"
			p.contexts[i.noncg[bases2 == "A" & bases3 != "G"]] <- "CAH"
			p.contexts[i.noncg[bases2 == "T" & bases3 == "G"]] <- "CTG"
			p.contexts[i.noncg[bases2 == "T" & bases3 != "G"]] <- "CTH"
		}
		result[["Context"]] <- p.contexts
		result
	}
#	l.neighborhoods <- as.list(rep(LENGTH.NEIGHBORHOOD, length(locations)))
#	result2 <- foreach(chrom = as.list(names(locations)), loci = locations, l.neighborhood = l.neighborhoods,
#		.packages = c("IRanges", "Biostrings", "GenomicRanges"), export = c('genome.data')) %dopar%
#		calculate.cg(chrom, loci, l.neighborhood)
#	names(result) <- names(locations)
	result <- list()
	for (chrom in names(locations)) {
		result[[chrom]] <- calculate.cg(chrom, locations[[chrom]], LENGTH.NEIGHBORHOOD)
	}

	## Copy the calculated values to probe.infos
	probe.infos[["CpG"]] <- as.integer(NA)
	probe.infos[["GC"]] <- as.integer(NA)
	probe.infos[["Context"]] <- factor("Other", levels = c("CG", "CC", "CAG", "CAH", "CTG", "CTH", "Other"))
	for (chromosome in names(result)) {
		i <- which(probe.infos[, "Chromosome"] == chromosome)
		for (cname in colnames(result[[chromosome]])) {
			probe.infos[i, cname] <- result[[chromosome]][[cname]]
		}
	}
	probe.infos[grep("^rs", probe.infos[["ID"]]), "Context"] <- "Other"

	logger.status("Calculated CpG density, GC content and probe target context")
	return(probe.infos)
}

########################################################################################################################

#' rnb.update.probe.annotation.msnps
#'
#' ...
#'
#' @param probe.infos ...
#' @param snps        ...
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.msnps <- function(probe.infos, snps = .globals[['snps']]) {
	alleles.A <- as.character(probe.infos[, "AlleleA Probe Sequence"])
	alleles.B <- as.character(probe.infos[, "AlleleB Probe Sequence"])
	bisulfite.convert <- function(x) {
		x.chars <- strsplit(x, split = NULL)
		sapply(x.chars, function(x.char) {
			for (j in 1:length(x.char)) {
				if (x.char[j] == "C" && (j == length(x.char) || x.char[j + 1] != "G")) {
					x.char[j] <- "T"
				}
			}
			paste(x.char, collapse = "")
		})
	}
	REV.COMPLEMENT <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
	reverse.complement <- function(x) {
		x.chars <- strsplit(x, split = NULL)
		sapply(x.chars, function(x.char) { paste(rev(REV.COMPLEMENT[x.char]), collapse = "") })
	}
	mismatches <- function(expected, observed) {
		mapply(function(x, y) sum(x != y), strsplit(expected, ""), strsplit(observed, ""))
	}
	genome.data <- get.genome.data()
	probe.regs <- list()
	for (chromosome in .globals[['CHROMOSOMES']]) {
#chromosome <- .globals[['CHROMOSOMES']][1]
		inds <- which((probe.infos[["Chromosome"]] == chromosome) & (probe.infos[["Strand"]] != "*"))
		loci <- probe.infos[inds, "Location"]
		probe.regs[[chromosome]] <- data.frame(
			"target" = loci, "start" = as.integer(NA), "end" = as.integer(NA), "strand" = probe.infos[inds, "Strand"])
		rownames(probe.regs[[chromosome]]) <- inds
		chrom.sequence <- genome.data[[chromosome]]

		alleles.A.exp <- rep(as.character(NA), times = length(inds))
		alleles.B.exp <- rep(as.character(NA), times = length(inds))

		## For type I forward:
		i <- which(probe.infos[inds, "Design"] == "I" & probe.infos[inds, "Strand"] == "+")
		if (length(i) != 0) {
			probe.regs[[chromosome]][i, "start"] <- loci[i]
			probe.regs[[chromosome]][i, "end"] <- loci[i] + 50L
			dna.seq <- as.character(suppressWarnings(Views(chrom.sequence, start = loci[i], end = loci[i] + 49L)))
			alleles.exp <- reverse.complement(bisulfite.convert(dna.seq))
			alleles.A.exp[i] <- gsub("G", "A", alleles.exp, fixed = TRUE)
			alleles.B.exp[i] <- alleles.exp
		}
		## For type I reverse:
		i <- which(probe.infos[inds, "Design"] == "I" & probe.infos[inds, "Strand"] == "-")
		if (length(i) != 0) {
			probe.regs[[chromosome]][i, "start"] <- loci[i] - 48L
			probe.regs[[chromosome]][i, "end"] <- loci[i] + 2L
			dna.seq <- as.character(suppressWarnings(Views(chrom.sequence, start = loci[i] - 48L, end = loci[i] + 1L)))
			alleles.exp <- reverse.complement(bisulfite.convert(reverse.complement(dna.seq)))
			alleles.A.exp[i] <- gsub("G", "A", alleles.exp, fixed = TRUE)
			alleles.B.exp[i] <- alleles.exp
		}
		## For type II forward:
		i <- which(probe.infos[inds, "Design"] == "II" & probe.infos[inds, "Strand"] == "+")
		if (length(i) != 0) {
			probe.regs[[chromosome]][i, "start"] <- loci[i] + 1L
			probe.regs[[chromosome]][i, "end"] <- loci[i] + 51L
			dna.seq <- as.character(suppressWarnings(Views(chrom.sequence, start = loci[i] + 1L, end = loci[i] + 50L)))
			alleles.exp <- reverse.complement(bisulfite.convert(dna.seq))
			alleles.A.exp[i] <- gsub("G", "R", alleles.exp, fixed = TRUE)
		}
		## For type II reverse:
		i <- which(probe.infos[inds, "Design"] == "II" & probe.infos[inds, "Strand"] == "-")
		if (length(i) != 0) {
			probe.regs[[chromosome]][i, "start"] <- loci[i] - 49L
			probe.regs[[chromosome]][i, "end"] <- loci[i] + 1L
			dna.seq <- as.character(suppressWarnings(Views(chrom.sequence, start = loci[i] - 49L, end = loci[i])))
			alleles.exp <- reverse.complement(bisulfite.convert(reverse.complement(dna.seq)))
			alleles.A.exp[i] <- gsub("G", "R", alleles.exp, fixed = TRUE)
		}

		probe.infos[inds, "Mismatches A"] <- mismatches(alleles.A.exp, alleles.A[inds])
		probe.infos[inds, "Mismatches B"] <- mismatches(alleles.B.exp, alleles.B[inds])
	}
	rm(alleles.A, alleles.B, bisulfite.convert, REV.COMPLEMENT, reverse.complement, mismatches)
	rm(chromosome, inds, loci, chrom.sequence, alleles.A.exp, alleles.B.exp, i)
	suppressWarnings(rm(dna.seq, alleles.exp))
	logger.status("Counted mismatches between expected and defined allele sequence")

	snp.stats <- foreach(pr.regs = probe.regs, snp.regs = snps[names(probe.regs)], .combine = rbind,
		.export = "rnb.update.probe.annotation.snps") %dopar%
		rnb.update.probe.annotation.snps(pr.regs, snp.regs)
	i <- as.integer(rownames(snp.stats))
	for (cname in colnames(snp.stats)) {
		probe.infos[[cname]] <- as.integer(NA)
		probe.infos[i, cname] <- snp.stats[, cname]
	}
	logger.status("Marked overlaps of probes with dbSNP records")

	probe.infos
}

########################################################################################################################

#' rnb.update.probe.annotation.snps
#'
#' Computes, for one chromosome, overlap between Infinium probe target locations and sequences with SNPs.
#'
#' @param probe.regs Infinium probe locations in a \code{data.frame}.
#' @param snps       Information about SNPs in the form of a \code{data.frame} containing at least the following
#'                   columns: \code{"start"}, \code{"end"}, \code{"C2T"}, \code{"G2A"}.
#' @return Matrix of three columns - \code{"SNPs 3"}, \code{"SNPs 5"} and \code{"SNPs Full"} - containing information on
#'         overlap of target dinucleotides and probe sequences with SNPs.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.snps <- function(probe.regs, snps) {
	snp.stats <- matrix(0L, nrow = nrow(probe.regs), ncol = 3,
		dimnames = list(rownames(probe.regs), c("SNPs 3", "SNPs 5", "SNPs Full")))

	if (!is.null(snps)) {
		## Count SNPs in probe regions
		for (i in 1:nrow(probe.regs)) {
			i.snps <- which(probe.regs[i, "start"] <= snps[, "end"] & snps[, "start"] <= probe.regs[i, "end"])
			if (length(i.snps) != 0) {
				cname <- ifelse(probe.regs[i, "strand"] == "-", "G2A", "C2T")
				i.snps <- i.snps[snps[i.snps, cname] == FALSE]
				if (length(i.snps) != 0) {
					if (probe.regs[i, "strand"] == "-") {
						p.ends <- c(probe.regs[i, "end"], probe.regs[i, "end"])
						p.starts <- p.ends - c(2L, 4L)
					} else {
						p.starts <- c(probe.regs[i, "start"], probe.regs[i, "start"])
						p.ends <- p.starts + c(2L, 4L)
					}
					snps.local <- snps[i.snps, ]
					snp.stats[i, 1] <- sum(p.starts[1] <= snps.local[, "end"] & snps.local[, "start"] <= p.ends[1])
					snp.stats[i, 2] <- sum(p.starts[2] <= snps.local[, "end"] & snps.local[, "start"] <= p.ends[2])
					snp.stats[i, 3] <- length(i.snps)
				}
			}
		}
	}

	return(snp.stats)
}

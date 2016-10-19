########################################################################################################################
## probes.R
## created: 2015-04-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Common functions for creating and updating Infinium probe annotation tables.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.get.illumina.annotation.columns
#'
#' Gets the expected column names from an Illumina's probe annotation table.
#'
#' @param assay The methylation assay: one of \code{"27k"}, \code{"450k"} or \code{"EPIC"}.
#' @return Columns to expect in the probe annotation table in the form of a \code{character} \code{vector}.
#' @author Yassen Assenov
#' @noRd
rnb.get.illumina.annotation.columns <- function(assay) {
	tbl <- read.csv(system.file("extdata/probeAnnotationColumns.csv", package = "RnBeadsAnnotationCreator"),
		check.names = FALSE, stringsAsFactors = FALSE)
	tbl <- tbl[!is.na(tbl[, assay]), c(assay, "Name", "RnBeads")]
	tbl <- tbl[order(tbl[, 1]), ]
	result <- tbl[, "RnBeads"]
	names(result) <- tbl[, "Name"]
	result
}

########################################################################################################################

#' rnb.load.probe.annotation.geo
#'
#' Loads the specified probe annotation table from GEO.
#'
#' @param ftp.table     FTP link to the probe annotation table in GEO.
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @param platform      Assay name; must be one \code{"HumanMehylation27"} or \code{"HumanMehylation450"}.
#' @return \code{list} of two \code{data.frame}s:
#'         \describe{
#'           \item{\code{"probes"}}{Probe annotation; column names are the ones specified in \code{table.columns}.}
#'           \item{\code{"controls"}}{Control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.load.probe.annotation.geo <- function(ftp.table, table.columns, platform) {

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
#' @param platform Assay name, must be one of \code{"HumanMethylation27"} or \code{"HumanMethylation450"}.
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
rnb.update.probe.annotation.methylumi <- function(platform = "HumanMethylation27") {

	pltrf <- ifelse(platform == "HumanMethylation27", "HM27", "HM450")
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
	if (platform == "HumanMethylation27") {
		control.probe.infos <- hm27.controls # [, "Address", "Type", "Color_Channel", "Name"]
	} else { # platform == "HumanMethylation450"
		control.probe.infos <- hm450.controls # [, "Address", "Type", "Color_Channel", "Name"]
	}
	rownames(control.probe.infos) <- control.probe.infos[, "Address"]
	logger.status("Extracted control probe annotation from the methylumi package")

	return(list("probes" = probe.infos, "controls" = control.probe.infos))
}

########################################################################################################################

#' rnb.probes.fix.infinium.columns
#' 
#' Validates the column values in an Inifinium probe annotation table.
#' @param probe.infos Probe annotation table to be validated.
#' @return The modified and sorted annotation table.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.probes.fix.infinium.columns <- function(probe.infos) {
	probe.infos[, "Chromosome"] <- as.character(probe.infos[, "Chromosome"])
	i <- which(probe.infos[, "Chromosome"] == "")
	if (length(i) != 0) { probe.infos[i, "Chromosome"] <- NA }
	i <- is.na(probe.infos[, "Chromosome"])
	if (!identical(i, is.na(probe.infos[["Location"]]))) {
		logger.error("Inconsistent values in columns Chromosome and Location")
	}
	if ("Design" %in% colnames(probe.infos)) {
		if (any(probe.infos[!i, "Strand"] == "" | probe.infos[!i, "Design"] == "")) {
			logger.error("Missing design and/or strand information for some probes with specified location")
		}
	}
	i <- which(!(i | grepl("^chr", probe.infos[, "Chromosome"])))
	if (length(i) != 0) {
		probe.infos[i, "Chromosome"] <- paste0("chr", probe.infos[i, "Chromosome"])
	}
	if (!setequal(names(.globals[['CHROMOSOMES']]), unique(na.omit(probe.infos[, "Chromosome"])))) {
		logger.error("Unexpected chromosome names in the probe definition table")
	}
	probe.infos[, "Chromosome"] <- factor(as.character(probe.infos[, "Chromosome"]),
			levels = names(.globals[['CHROMOSOMES']]))
	if ("Color" %in% colnames(probe.infos)) {
		if (!identical(levels(probe.infos[, "Color"]), c("", "Grn", "Red"))) {
			logger.error("Unexpected color channel values in the probe definition table")
		}
		levels(probe.infos[, "Color"]) <- c("Both", "Grn", "Red")
	}
	probe.infos[["ID"]] <- as.character(probe.infos[["ID"]])
	rownames(probe.infos) <- probe.infos[["ID"]]
	if ("Random" %in% colnames(probe.infos)) {
		probe.infos[is.na(probe.infos[["Random"]]), "Random"] <- FALSE
	}
	if ("HumanMethylation27" %in% colnames(probe.infos)) {
		probe.infos[is.na(probe.infos[["HumanMethylation27"]]), "HumanMethylation27"] <- FALSE
	}
	if ("Strand" %in% colnames(probe.infos)) {
		probe.infos[, "Strand"] <- rnb.fix.strand(probe.infos[, "Strand"])
	}

	## Improve the notation of CGI relation
	cgi.relations <- c("Open Sea", "Island", "North Shelf", "North Shore", "South Shelf", "South Shore")
	names(cgi.relations) <- c("", "Island", "N_Shelf", "N_Shore", "S_Shelf", "S_Shore")
	if (!identical(levels(probe.infos[, "CGI Relation"]), names(cgi.relations))) {
		logger.error("Unexpected values in column for relation to CpG island")
	}
	levels(probe.infos[, "CGI Relation"]) <- cgi.relations

	## Sort based on chromosome and position
	probe.infos[with(probe.infos, order(Chromosome, Location)), ]
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

	genome.data <- rnb.genome.data()
	calculate.cg <- function(chrom, loci, l.neighborhood) {
		starts <- loci - l.neighborhood / 2L + 1L
		ends <- loci + l.neighborhood / 2L
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
	genome.data <- rnb.genome.data()
	probe.regs <- list()
	for (chromosome in names(.globals[['CHROMOSOMES']])) {
#chromosome <- names(.globals[['CHROMOSOMES']])[1]
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

	if (!is.null(snps)) {
		snp.stats <- suppressWarnings(
			foreach(pr.regs = probe.regs, snp.regs = snps[names(probe.regs)], .combine = rbind,
					.export = "rnb.update.probe.annotation.snps.chrom") %dopar%
				rnb.update.probe.annotation.snps.chrom(pr.regs, snp.regs))
		i <- as.integer(rownames(snp.stats))
		for (cname in colnames(snp.stats)) {
			probe.infos[[cname]] <- as.integer(NA)
			probe.infos[i, cname] <- snp.stats[, cname]
		}
		logger.status("Marked overlaps of probes with dbSNP records")
	}

	probe.infos
}

########################################################################################################################

#' rnb.update.probe.annotation.snps
#'
#' Computes the overlaps between Infinium probe target locations and sequences with SNPs.
#'
#' @param probe.infos Infinium probe annotation table as a \code{data.frame} containig at least the following columns:
#'                    \code{"Chromosome"}, \code{"Location"}, \code{"Design"}, \code{"Strand"}.
#' @param snps        Information about SNPs in the form of a \code{list} of \code{data.frame}s, each containing at
#'                    least the following columns: \code{"start"}, \code{"end"}, \code{"C2T"}, \code{"G2A"}.
#' @return The updated probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.snps <- function(probe.infos, snps = .globals[['snps']]) {
	if (is.null(snps)) {
		return(probe.infos)
	}

	probe.regs <- data.frame(
		"start" = as.integer(NA), "end" = as.integer(NA),
		"strand" = probe.infos[, "Strand"], "target" = probe.infos[, "Location"])
	genome.data <- rnb.genome.data()
	for (chromosome in names(.globals[['CHROMOSOMES']])) {
		chrom.sequence <- genome.data[[chromosome]]
		for (pr.design in c("I", "II")) {
			for (pr.strand in c("+", "-")) {
				i <- which(probe.infos[["Chromosome"]] == chromosome & probe.infos[["Design"]] == pr.design &
						 probe.infos[["Strand"]] == pr.strand)
				if (length(i) != 0) {
					p.coords <- rnb.infinium.probe.coords(probe.infos[i, "Location"], pr.design, pr.strand)
					probe.regs[i, "start"] <- p.coords[, 1]
					probe.regs[i, "end"] <- p.coords[, 2]
					rm(p.coords)
				}
			}
		}
	}
	rm(genome.data, chromosome, chrom.sequence, pr.design, pr.strand, i)

	probe.regs <- tapply(1:nrow(probe.regs), probe.infos$Chromosome, function(i) { probe.regs[i, ] })
	snp.stats <- suppressWarnings(
		foreach(pr.regs = probe.regs, snp.regs = snps[names(probe.regs)], .combine = rbind,
				.export = "rnb.update.probe.annotation.snps.chrom") %dopar%
			rnb.update.probe.annotation.snps.chrom(pr.regs, snp.regs))
	i <- as.integer(rownames(snp.stats))
	for (cname in colnames(snp.stats)) {
		probe.infos[[cname]] <- as.integer(NA)
		probe.infos[i, cname] <- snp.stats[, cname]
	}
	logger.status("Marked overlaps of probes with dbSNP records")
	probe.infos
}

########################################################################################################################

#' rnb.update.probe.annotation.snps.chrom
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
rnb.update.probe.annotation.snps.chrom <- function(probe.regs, snps) {
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

########################################################################################################################

#' rnb.update.probe.annotation.cr
#'
#' Extracts information on cross-reactive probes.
#'
#' @param probe.ids All probe identifiers in a \code{character} vector.
#' @param platform      Assay name; must be one \code{"HumanMehylation27"} or \code{"HumanMehylation450"}.
#' @return \code{integer} vector of the same length as \code{probe.ids}. The \code{i}th element in this vector stores
#'         the number of cross-reactive targets matching 47 or more bases in the sequence of probe \code{i}.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.cr <- function(probe.ids, platform) {
	if (platform == "HumanMethylation27") {
		fname <- "probes27"
	} else if (platform == "HumanMethylation450") {
		fname <- "probes450"
	} else { # platform == "MethylationEPIC"
		fname <- "probesEPIC"
	}
	fname <- paste0("extdata/", .globals[['assembly']], ".", fname, ".crossreactive.txt")
	fname <- system.file(fname, package = "RnBeadsAnnotationCreator")
	if (file.exists(fname)) {
		tbl <- read.delim(fname, quote = "", check.names = FALSE, stringsAsFactors = FALSE)
		logger.status(c("Loaded list of cross-reactive probes from", fname))
		if (nrow(tbl) < 1) {
			logger.error("Empty table")
		}
		expected <- c("TargetID" = "character", "47" = "integer", "48" = "integer", "49" = "integer", "50" = "integer")
		if (!identical(sapply(tbl, class), expected)) {
			logger.error("Unexpected table structure")
		}
		if (anyDuplicated(tbl[, 1]) != 0) {
			logger.error("Unexpected table structure; duplicated probe identifiers found")
		}
		rownames(tbl) <- tbl[, 1]
		tbl <- as.matrix(tbl[, -1])
		if (any(tbl < 0L)) {
			logger.error("Unexpected table structure; negative counts found")
		}
		tbl <- apply(tbl, 1, sum, na.rm = TRUE)
		ids.notfound <- setdiff(names(tbl), probe.ids)
		if (length(ids.notfound) != 0) {
			if (length(ids.notfound) == length(tbl)) {
				logger.error("None of the probe IDs in the table is a valid probe ID")
			}
			logger.warning(c("Ignoring", length(ids.notfound), "invalid probe IDs in the table"))
			tbl <- tbl[intersect(probe.ids, names(tbl))]
		}
		result <- rep(0L, length(probe.ids))
		names(result) <- probe.ids
		result[names(tbl)] <- tbl
		result <- unname(result)
		logger.status("Added information on cross-reactive probes")
	} else {
		result <- rep(as.integer(NA), length(probe.ids))
		logger.warning("No list of cross-reactive probes available; creating a column with NAs")
	}
	result
}

########################################################################################################################

#' rnb.probe.infos.to.GRanges
#'
#' Converts a data frame with probe annotation to a \code{GRangesList} instance.
#'
#' @param probe.infos Probe annotation as a \code{data.frame} containing at least the following columns:
#'                    \code{"Chromosome"}, \code{"Location"}, \code{"ID"}.
#' @return Converted annotation in the form of a \code{GRangesList} object, one \code{GRanges} instance per chromosome.
#' @author Yassen Assenov
#' @noRd
rnb.probe.infos.to.GRanges <- function(probe.infos) {
	starts <- probe.infos[, "Location"]
	starts[is.na(starts)] <- 0L
	cnames <- c("Strand", "AddressA", "AddressB", "Design", "Color", "Context", "Random", "HumanMethylation27",
		"HumanMethylation450", "Mismatches A", "Mismatches B", "CGI Relation", "CpG", "GC", "SNPs 3", "SNPs 5",
		"SNPs Full", "Cross-reactive")
	cnames <- as.list(probe.infos[, intersect(cnames, colnames(probe.infos))])
	result <- c(list(seqnames = probe.infos[, "Chromosome"],
			ranges = IRanges(start = starts, end = starts + 1L, names = probe.infos[, "ID"]),
			strands <- rnb.fix.strand(probe.infos[, "Strand"])),
		cnames, list(seqinfo = seqinfo(rnb.genome.data())[names(.globals[['CHROMOSOMES']]), ]))
	result <- rnb.sort.regions(do.call(GRanges, result))
	result <- GenomicRanges::split(result, seqnames(result))[names(.globals[['CHROMOSOMES']])]
	seqinfo(result) <- seqinfo(rnb.genome.data())[names(.globals[['CHROMOSOMES']]), ]
	result
}

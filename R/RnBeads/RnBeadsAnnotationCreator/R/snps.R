########################################################################################################################
## snps.R
## created: 2012-08-18
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes tables of SNP annotations and enriches the probe annotation tables with information about overlapping
## SNPs.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.download.dbsnp
#'
#' Downloads all required VCF files from the latest release of dbSNP.
#'
#' @param ftp.files Full FTP paths of the VCF files in dbSNP.
#' @param base.dir  Local directory to store the downloaded VCF files.
#' @return \code{character} vector of length \code{length(ftp.files)}, storing the names of all local copies of the
#'         downloaded VCF files.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.download.dbsnp <- function(ftp.files, base.dir) {
	if (!(is.character(ftp.files) && length(ftp.files) > 0 && all(!is.na(ftp.files)))) {
		stop("invalid value for ftp.files; expected character")
	}
	base.dir <- base.dir[1]

	fnames <- gsub("^.+/([^/]+)$", "\\1", ftp.files)
	dest.files <- file.path(base.dir, fnames)
	for (i in 1:length(ftp.files)) {
		if (file.exists(dest.files[i])) {
			logger.info(c("File", fnames[i], "already downloaded"))
		} else {
			logger.start(c("Downloading", fnames[i]))
			if (download.file(ftp.files[i], dest.files[i], quiet = TRUE, mode = "wb") != 0) {
				logger.error(c("Could not download", fnames[i]))
			}
			logger.completed()
		}
	}
	dest.files
}

########################################################################################################################

#' rnb.update.load.vcf
#'
#' Loads information on the variations from the given VCF file. downloaded from dbSNP.
#'
#' @param fname Name of the VCF file to load.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.load.vcf <- function(fname) {

	## Extract meta information and header
	txt <- scan(fname, "", sep = "\n", quiet = TRUE)
	logger.status(c("Loaded", length(txt), "lines from", fname))
	i.header <- grep("^#", txt)
	meta.regex <- "^##([^=]+)=(.+)$"
	i.meta <- grep(meta.regex, txt[i.header])
	if (!(length(i.meta) > 0 && identical(i.header, 1:length(i.header)) && identical(i.meta, 1:length(i.meta)))) {
		logger.error("Unsupported structure of the VCF file")
	}
	meta.names <- gsub(meta.regex, "\\1", txt[i.meta])
	meta.values <- gsub(meta.regex, "\\2", txt[i.meta])

	## Validate header information
	g.value <- function(m.name, required = TRUE) {
		i <- which(meta.names == m.name)
		if (length(i) == 0) {
			if (required) {

			}
			return(NULL)
		}
		if (length(i) != 1) {
			msg <- ifelse(length(i) == 0, "missing", "multiple values in")
			stop(msg, " meta information for ", m.name)
		}
		meta.values[i]
	}
	if (g.value("fileformat") != "VCFv4.0") {
		logger.error("unsupported VCF format, expected VCFv4.0")
	}
	ref <- g.value("reference")
	if (!(ref %in% names(REFERENCE2ASSEMBLY))) {
		logger.error(c("Unsupported reference genome:", ref))
	}
	if (REFERENCE2ASSEMBLY[ref] != .globals[['assembly']]) {
		logger.error(c("Invalid genome assembly:", REFERENCE2ASSEMBLY[ref], ", expected", .globals[['assembly']]))
	}
	version.string <- g.value("fileDate")
	if (is.null(version.string)) {
		version.string <- ""
	} else {
		version.string <- gsub("^(\\d{4})(\\d{2})(\\d{2})$", " from \\1-\\2-\\3", version.string)
	}
	version.string <- paste0("dbSNP ", g.value("dbSNP_BUILD_ID"), version.string, ", reference ", ref)
	cnames <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
	if (!identical(strsplit(txt[length(i.header)], "\t", fixed = TRUE)[[1]], cnames)) {
		logger.error("Unexpected columns in the VCF table")
	}

	## Extract SNP information
	txt <- strsplit(txt[-(1:length(i.header))], "\t", fixed = TRUE)
	i <- which(sapply(txt, length) != length(cnames))
	if (length(i) != 0) {
		logger.error(c("unexpected number of columns at line", (i[1] + length(i.header))))
	}
	# remove duplicate rows
	dupRows <- duplicated(txt)
	if(sum(dupRows) > 0){
		txt <- txt[!dupRows]
		logger.warning(c("Removed",sum(dupRows),"duplicate rows"))
	}
	chroms <- sapply(txt, '[', 1)
	CHROMOSOMES <- .globals[['CHROMOSOMES']]
	is.valid.chromosome <- (chroms %in% names(CHROMOSOMES))
	if (all(!is.valid.chromosome)) {
		# FIXME: Use a general mapping from short names to long names in chromosomes; "MT" should map to "chrM"
		i <- which(chroms %in% unname(CHROMOSOMES))
		if (length(i) != 0) {
			short2long.name <- names(CHROMOSOMES)
			names(short2long.name) <- unname(CHROMOSOMES)
			chroms[i] <- unname(short2long.name[chroms[i]])
			rm(short2long.name)
		}
		is.valid.chromosome <- (chroms %in% names(CHROMOSOMES))
	}
	logger.info(c("Records with supported chromosome:", sum(is.valid.chromosome), ", with unsupported ones:",
		sum(!is.valid.chromosome)))
	ids <- sapply(txt, '[', 3)
	i <- anyDuplicated(ids)
	if (i != 0) {
		logger.error(paste0("duplicated identifier (", ids[i], ") found at line ", (i[1] + length(i.header))))
	}
	pos <- suppressWarnings(as.integer(sapply(txt, '[', 2)))
	i <- which(is.na(pos))
	if (length(i) != 0) {
		logger.error(c("invalid genomic position at line", (i[1] + length(i.header))))
	}
	infos <- sapply(txt, '[', 8)

	## Extract allele origin
	regex.ao <- "^.*SAO=([0-3]).*$" # 0 - unspecified, 1 - germline, 2 - somatic, 3 - both
	i <- which(!grepl(regex.ao, infos))
	if (length(i) != 0) {
		logger.error(c("missing or invalid variant allele origin at line", (i[1] + length(i.header))))
	}
	allele.origin <- as.integer(gsub(regex.ao, "\\1", infos))
	is.valid.allele <- allele.origin != 2L
	logger.info(c("Records with supported allele origin:", sum(is.valid.allele), ", with unsupported ones:",
		sum(!is.valid.allele)))

	## Extract major allele frequencies
	regex.frequency <- "^.*CAF=\\[([0-9\\.,]+)\\].*$"
	i <- which(grepl(regex.frequency, infos))
	if (length(i) == 0) {
		major.frequency <- rep(as.double(NA), length(infos))
		logger.warning(c("Could not detect MAF data (CAF field)"))
		regex.frequency.g5ind <- "(^G5.*$|^.*;G5A?$|^.*;G5A?;.*$)"
		hasG5tag <- grepl(regex.frequency.g5ind, infos)
		if (sum(hasG5tag) > 0) {
			is.valid.frequency <- hasG5tag
		} else {
			is.valid.frequency <- rep(TRUE, length(infos))
			logger.warning(c("Could not detect MAF data (G5 field); allele frequency is not considered"))
		}
	} else {
		major.frequency <- rep(as.double(NA), length(infos))
		frequencies <- strsplit(gsub(regex.frequency, "\\1", infos[i]), ",", fixed = TRUE)
		major.frequency[i] <- suppressWarnings(sapply(frequencies, function(x) { max(as.double(x), na.rm = TRUE) }))
		rm(frequencies)
		is.valid.frequency <- (major.frequency <= MAJOR.ALLELE.FREQUENCY)
		logger.info(c("Records with supported MAF:", sum(is.valid.frequency, na.rm = TRUE), ", with unsupported ones:",
			sum(!is.valid.frequency, na.rm = TRUE), ", with unknown:", sum(is.na(major.frequency))))
		is.valid.frequency[is.na(is.valid.frequency)] <- FALSE
	}

	## Filter based on chromosome name, allele origin and major allele frequency
	i <- which(is.valid.chromosome & is.valid.allele & is.valid.frequency)
	result <- data.frame(chromosome = chroms[i], location = pos[i],
		base.reference = sapply(txt, '[', 4)[i], base.alt = sapply(txt, '[', 5)[i],
		origin = allele.origin[i], frequency = major.frequency[i],
		row.names = ids[i], check.names = FALSE, stringsAsFactors = FALSE)
	result$chromosome <- factor(result$chromosome, levels = names(CHROMOSOMES))
	result <- result[with(result, order(chromosome, location)), ]
	attr(result, "version") <- version.string
	result
}

########################################################################################################################

#' rnb.construct.snp.types
#'
#' Processes the downloaded dbSNP records and constructs a table with polymorphism types.
#'
#' @param dnp.df \code{data.frame} with downloaded dbSNP records for one chromosome. This function uses the columns
#'               \code{"location"}, \code{"base.reference"} and \code{"base.alt"}.
#' @return \code{data.frame} with polymorphism records sorted by starting position.
#'
#' @author Yassen Assenov
#' @noRd
rnb.construct.snp.types <- function(snp.df) {
	result <- data.frame(
		start = rep(0L, nrow(snp.df)),
		end = rep(0L, nrow(snp.df)),
		type = factor("replacement", levels = c("deletion", "deletion+replacement", "replacement", "insertion", "other")),
		C2T = FALSE,
		G2A = FALSE,
		row.names = rownames(snp.df))
	alt.sequences <- strsplit(snp.df[, "base.alt"], ",", fixed = TRUE)
	for (i in 1:nrow(snp.df)) {
		i.start <- snp.df[i, "location"]
		i.sequences <- c(snp.df[i, "base.reference"], alt.sequences[[i]])
		i.nchar <- nchar(i.sequences)
		while (all(i.nchar != 0)) {
			if (length(unique(substr(i.sequences, 1L, 1L))) == 1L) {
				i.sequences <- substring(i.sequences, 2L)
				i.start <- i.start + 1L
				i.nchar <- i.nchar - 1L
			} else {
				break
			}
		}
		result[i, 1:2] <- c(i.start, i.start + i.nchar[1] - 1L)
		i.nchar <- c(i.nchar[1], range(i.nchar[-1]))
		if (i.nchar[1] > i.nchar[3]) {
			if (i.nchar[3] == 0) {
				result[i, 3] <- "deletion"
			} else {
				result[i, 3] <- "deletion+replacement"
			}
		} else if (i.nchar[1] < i.nchar[2]) {
			if (i.nchar[1] == 0) {
				result[i, 1:2] <- result[i, 2:1]
				result[i, 3] <- "insertion"
			} else {
				## Record combines multiple modifications; set deletion+replacement as a workaround
				result[i, 3] <- "deletion+replacement"
			}
		} else if (i.nchar[1] == i.nchar[2] && i.nchar[2] == i.nchar[3]) {
			if (identical("C", i.sequences[1]) && identical("T", i.sequences[-1])) {
				result[i, "C2T"] <- TRUE
			} else if (identical("G", i.sequences[1]) && identical("A", i.sequences[-1])) {
				result[i, "G2A"] <- TRUE
			}
		} else {
			result[i, 3] <- "other"
		}
	}
	result[order(result[, 1], result[, 2]), ]
}

########################################################################################################################

#' rnb.update.dbsnp
#'
#' Downloads and processes the information from dbSNP.
#'
#' @param ftp.files ...
#' @return SNP annotation in the form of a \code{list} of \code{data.frame} objects, one per chromosome.
#' @author Yassen Assenov
#' @noRd
rnb.update.dbsnp <- function(ftp.files) {
	fname <- file.path(.globals[["DIR.PACKAGE"]], "temp", "snps.all.RData")
	if (file.exists(fname)) {
		load(fname) # -> snps, db.version
		logger.status(c("Loaded SNP tables from", fname))
	} else {
		base.dir <- file.path(.globals[["DIR.PACKAGE"]], "temp", "snps")
		if (!file.exists(base.dir)) {
			if (!dir.create(base.dir, recursive = TRUE)) {
				logger.error(c("Could not create directory", base.dir))
			}
		} else if (!isTRUE(file.info(base.dir)[, "isdir"])) {
			logger.error(c("Expected directory", base.dir, "is a regular file"))
		}
		fnames <- rnb.update.download.dbsnp(ftp.files, base.dir)
		logger.start(paste0("Loading Downloaded File", ifelse(length(fnames) != 1, "s", "")))
		snps <- lapply(fnames, rnb.update.load.vcf)
		logger.completed()

		db.version <- sapply(snps, attr, "version")
		if (length(unique(db.version)) != 1) {
			logger.warning(c("Differing version strings in the downloaded files; versions are concatenated"))
		}
		db.version <- paste(unique(db.version), collapse = " ; ")
		snps <- do.call(rbind, snps)
		if (all(is.na(snps$frequency))) {
			snps$frequency <- NULL
		}
		snps <- snps[with(snps, order(chromosome, location)), ]
		snps <- tapply(1:nrow(snps), snps$chromosome, function(i) { snps[i, ] })
		logger.info(c("Total number of SNP records:", sum(sapply(snps, nrow))))

		save(snps, db.version, file = fname, compression_level = 9L)
		logger.status(c("Saved SNP tables to", fname))
	}

	## Construct tables of polymorphism types
	if (.globals[['SNP_MAX']] != 0 && .globals[['SNP_MAX']] <= max(sapply(snps, nrow))) {
		fname <- normalizePath(fname, '/')
		dname <- normalizePath(file.path(.globals[["DIR.PACKAGE"]], "temp", "qsub"), '/', FALSE)
		txt <- c('suppressPackageStartupMessages(library(RnBeadsAnnotationCreator))',
			paste0('txt <- rnb.split.snps("', fname, '", "', dname, '")'),
			'cat(txt, sep = "\\n", file = "split.snps.log")')
		fname <- file.path(.globals[["DIR.PACKAGE"]], "temp", "split.snps.R")
		cat(paste(txt, collapse = "\n"), file = fname)
		logger.error(paste('dbSNP table(s) are too large. MapReduce implemented in', fname))
	}

#	chromosomes <- names(snps)
#	snps <- suppressWarnings(foreach(snp.df = snps) %dopar% rnb.construct.snp.types(snp.df))
#	names(snps) <- chromosomes
	fnames <- file.path(.globals[["DIR.PACKAGE"]], "temp", "snps", paste0("snps.", names(snps), ".RDS"))
	names(fnames) <- names(snps)
	for (chrom in names(snps)) {
		fname <- fnames[chrom]
		if (file.exists(fname)) {
			snps[[chrom]] <- readRDS(fname)
			logger.status(c("Loaded", chrom))
		} else {
			snps[[chrom]] <- rnb.construct.snp.types(snps[[chrom]])
			con <- gzfile(fname, "wb", compression = 9L); saveRDS(snps[[chrom]], con); close(con); rm(con)
			logger.status(c("Processed", chrom))
		}
		invisible(gc())
	}
	attr(snps, "version") <- db.version
	logger.status("Constructed tables of polymorphism types")

	return(snps)
}

########################################################################################################################

#' rnb.split.snps
#'
#' Performs dbSNP record post-processing using the MapReduce programming model. This function defines and submits jobs
#' to a computational cluster using the PBS scheduling system.
#'
#' @param snps.file      Full path of the file (named \code{snps.RData}) that stores the combined and filtered table(s)
#'                       from dbSNP. The processed tables (one per chromosome) that result from calling this function
#'                       will be written to dedicated files in a subdirectory named \code{snps} of the directory in
#'                       which this file is located.
#' @param temp.directory Full path of the directory to be for the storage of intermediate data files, log files and
#'                       scripts. This must be a non-existent path as this function attempts to create it. Note that
#'                       this temporary directory is \emph{not} cleaned after the completion of all jobs.
#' @param R.executable   Full path of the R executable file; used for the execution in R scripts in the PBS system.
#' @param batch.size     Maximum size, in number of records, of a table processed in a single job. This must be an
#'                       \code{integer} value between \code{10^3} and \code{10^6}.
#' @return Invisibly, the identifiers of all processes submitted to PBS.
#'
#' @author Yassen Assenov
#' @export
rnb.split.snps <- function(snps.file, temp.directory, R.executable = paste0(Sys.getenv("R_HOME"), "/bin/R"),
	batch.size = 200000L) {

	## Validate parameters
	validate.path <- function(x, pname, is.file = TRUE) {
		if (!(is.character(x) && length(x) == 1 && isTRUE(x) != "")) {
			stop(paste("invalid value for", pname))
		}
		if (!grepl("^(/|[A-Za-z]:).+$", x)) {
			stop(paste0("invalid value for ", pname, "; absolute path expected"))
		}
		if (grepl('"', x, fixed = TRUE)) {
			stop(paste0("invalid value for ", pname, "; double quotes are not allowed"))
		}
		pexists <- file.exists(x)
		if (is.file) {
			if (!pexists) {
				stop(paste0("invalid value for ", pname, "; file not found"))
			}
		} else if (pexists) {
			stop(paste0("invalid value for ", pname, "; must be a non-existent path"))
		}
	}
	validate.path(snps.file, "snps.file")
	validate.path(temp.directory, "temp.directory", FALSE)
	if (is.double(batch.size) && isTRUE(all(batch.size == as.integer(batch.size)))) {
		batch.size <- as.integer(batch.size)
	}
	if (!(is.integer(batch.size) && length(batch.size) == 1)) {
		stop("invalid value for batch.size")
	}
	if (!isTRUE(1000L < batch.size && batch.size <= 1000000L)) {
		stop("invalid value for batch.size; expected a value between 10^3 and 10^6")
	}
	validate.path(R.executable, "R.executable")

	## Load the tables of SNP records
	result <- tryCatch(suppressWarnings(load(snps.file)), error = function(err) { character() })
	if (!setequal(result, c("snps", "db.version"))) {
		stop(paste("File", snps.file, "is inaccessible or invalid"))
	}
	if (!(is.list(snps) && length(snps) != 0 && all(sapply(snps, class) == "data.frame") && (!is.null(names(snps))))) {
		stop("Loaded object snps contains invalid data")
	}
	if (!isTRUE(all(grepl("^[A-Za-z0-9_]+$", names(snps))))) {
		stop("Loaded object snps contains invalid names")
	}
	rm(validate.path, result, db.version)

	## Attempt to create the temporary directory
	if (!dir.create(temp.directory, FALSE, TRUE)) {
		stop("could not create temp.directory")
	}
	output.directory <- file.path(dirname(snps.file), "snps")
	if (!file.exists(output.directory)) {
		if (!dir.create(output.directory, FALSE)) {
			unlink(temp.directory)
			stop(c("could not create", output.directory))
		}
	}

	## Split the SNP records into smaller tables
	batch.chromosomes <- rep(names(snps), as.integer(ceiling(sapply(snps, nrow) / batch.size)))
	batch.chromosomes <- factor(batch.chromosomes, levels = names(snps))
	i <- 0L
	for (chrom in names(snps)) {
		n <- 0L
		while (n != nrow(snps[[chrom]])) {
			n.start <- n + 1L
			n <- min(n + batch.size, nrow(snps[[chrom]]))
			i <- i + 1L
			fname <- sprintf("db.%05d.RDS", i)
			con <- gzfile(file.path(temp.directory, fname), "wb", compression = 9L)
			saveRDS(snps[[chrom]][n.start:n, ], con)
			close(con)
		}
	}
	rm(snps, i, chrom, n, n.start, fname, con)
	invisible(gc())

	## Generate an R script for processing a single file
	txt <- capture.output(show(RnBeadsAnnotationCreator:::rnb.construct.snp.types))
	txt[1] <- paste("rnb.construct.snp.types <-", txt[1])
	txt[length(txt)] <- ''
	txt <- c(txt, 'i <- commandArgs()[length(commandArgs())]',
		'fnames <- paste0(c("db.", "tp."), i, ".RDS")',
		'tbl <- rnb.construct.snp.types(readRDS(fnames[1]))',
		'con <- gzfile(fnames[2], "wb", compression = 9L)',
		'saveRDS(tbl, con)',
		'close(con)',
		'file.remove(fnames[1])', '')
	cat(paste(txt, collapse = "\n"), file = file.path(temp.directory, "snp.types.R"))

	## Generate a shell script for processing a single file
	generate.shell <- function(fname, logfile, R.args, mem = NULL, walltime = NULL) {
		txt <- c('#!/bin/bash', '')
		if (!is.null(mem)) { txt <- c(txt, paste0('#PBS -l mem=', mem)) }
		if (!is.null(walltime)) { txt <- c(txt, paste0('#PBS -l walltime=', walltime, ':00:00')) }
		txt <- c(txt, '#PBS -l nodes=1:ppn=1',
			'#PBS -m ae',
			'#PBS -j oe',
			paste0('#PBS -o "', temp.directory, '/', logfile, '.log"'), '',
			paste0('cd "', temp.directory, '"'),
			paste0('"', R.executable, '" --no-restore --no-save --args ', R.args, ' < ', fname, '.R'), '')
		cat(paste(txt, collapse = "\n"), file = file.path(temp.directory, paste0(fname, ".sh")))
	}
	generate.shell('snp.types', '${i}', '${i}', paste0(ceiling(batch.size / 2000L), "m"),
		as.integer(ceiling(batch.size / 70000)))

	## Generate an R script for combining results
	txt <- c('arguments <- tail(commandArgs(), 3)',
		'ii <- as.integer(arguments[1:2])',
		'fnames <- sprintf("tp.%05d.RDS", ii[1]:ii[2])',
		'tbl <- do.call(rbind, unname(lapply(fnames, readRDS)))',
		'tbl <- tbl[order(tbl[, 1], tbl[, 2]), ]',
#		'dname <- dirname(arguments[3])',
#		'if (!file.exists(dname)) { dir.create(dname, FALSE, TRUE) }',
		'con <- gzfile(arguments[3], "wb", compression = 9L)',
		'saveRDS(tbl, con)',
		'close(con)',
		'fnames[!file.remove(fnames)]', '')
	fname <- file.path(temp.directory, "snp.combine.R")
	cat(paste(txt, collapse = "\n"), file = fname)

	## Generate a shell script for combining results
	generate.shell('snp.combine', '${chrom}', '${istart} ${iend} ${fn}')

	## Submit jobs for processing
	ids1 <- rep("", length(batch.chromosomes))
	for (i in 1:length(batch.chromosomes)) {
		ids1[i] <- sprintf('qsub -N SNPs_%05d -v i=%05d "%s/snp.types.sh"', i, i, temp.directory)
		ids1[i] <- system(ids1[i], intern = TRUE)
		Sys.sleep(1)
	}
	rm(generate.shell, i)

	## Submit jobs for combining
	ids2 <- character()
	for (chrom in levels(batch.chromosomes)) {
		i <- range(which(batch.chromosomes == chrom))
		fname <- paste0(output.directory, '/snps.', chrom, '.RDS')
		mem <- as.integer(ceiling(0.12 * (i[2] - i[1] + 1)))
		txt <- paste0('qsub -N SNPs_combine_', chrom, ' -l mem=', mem, 'g -l walltime=', mem, ':00:00 ',
			'-W depend=afterok:', paste0(ids1[i[1]:i[2]], collapse = ':'),
			' -v chrom=', chrom, ',istart=', i[1], ',iend=', i[2], ',fn="', fname, '" ',
			'"', temp.directory, '/snp.combine.sh"')
		ids2 <- c(ids2, system(txt, intern = TRUE))
		Sys.sleep(1)
	}
	invisible(c(ids1, ids2))
}

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

	reference2assembly <- c(
		"GRCh37.p10" = "hg19",
		"GCF_000001635.21" = "mm10", # Genome Reference Consortium Mouse Build 38 patch release 1 (GRCm38.p1) 
		"GCF_000001895.4" = "rn5")

	## Extract meta information and header
	txt <- scan(fname, "", sep = "\n", quiet = TRUE)
	logger.status(c("Loaded", length(txt), "lines from", fname))
	i.header <- grep("^#", txt)
	meta.regex <- "^##([^=]+)=(.+)$"
	i.meta <- grep(meta.regex, txt[i.header])
	if (!(length(i.meta) > 0 && identical(i.header, 1:length(i.header)) && identical(i.meta, 1:length(i.meta)))) {
		stop("unsupported structure of the VCF file")
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
		stop("unsupported VCF format, expected VCFv4.0")
	}
	ref <- g.value("reference")
	if (!(ref %in% names(REFERENCE2ASSEMBLY))) {
		stop(paste("unsupported reference genome:", ref))
	}
	if (REFERENCE2ASSEMBLY[ref] != assembly) {
		logger.error(c("Invalid genome assembly:", REFERENCE2ASSEMBLY[ref], ", expected", assembly))
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
		stop("unexpected columns in the VCF table")
	}

	## Extract SNP information
	txt <- strsplit(txt[-(1:length(i.header))], "\t", fixed = TRUE)
	i <- which(sapply(txt, length) != length(cnames))
	if (length(i) != 0) {
		stop("unexpected number of columns at line ", (i[1] + length(i.header)))
	}
	chroms <- sapply(txt, '[', 1)
	is.valid.chromosome <- (chroms %in% CHROMOSOMES)
	if (all(!is.valid.chromosome)) {
		chroms <- as.character(rnb.get.chromosomes(assembly)[chroms])
		is.valid.chromosome <- (chroms %in% CHROMOSOMES)
	}
	logger.info(c("Records with supported chromosome:", sum(is.valid.chromosome), ", with unsupported ones:",
		sum(!is.valid.chromosome)))
	ids <- sapply(txt, '[', 3)
	i <- anyDuplicated(ids) 
	if (i != 0) {
		stop("duplicated identifier (", ids[i], ") found at line ", (i[1] + length(i.header)))
	}
	pos <- suppressWarnings(as.integer(sapply(txt, '[', 2)))
	i <- which(is.na(pos))
	if (length(i) != 0) {
		stop("invalid genomic position at line ", (i[1] + length(i.header)))
	}
	ids <- sapply(txt, '[', 3)
	infos <- sapply(txt, '[', 8)
	
	## Extract allele origin
	regex.ao <- "^.*SAO=([0-3]).*$" # 0 - unspecified, 1 - Germline, 2 - Somatic, 3 - Both
	i <- which(!grepl(regex.ao, infos))
	if (length(i) != 0) {
		stop("missing or invalid variant allele origin at line ", (i[1] + length(i.header)))
	}
	allele.origin <- as.integer(gsub(regex.ao, "\\1", infos))
	is.valid.allele <- allele.origin != 2L
	logger.info(c("Records with supported allele origin:", sum(is.valid.allele), ", with unsupported ones:",
		sum(!is.valid.allele)))
	
	## Extract major allele frequencies
	regex.frequency <- "^.*CAF=\\[([0-9\\.,]+)\\].*$"
	i <- which(grepl(regex.frequency, infos))
	if (length(i) == 0) {
		is.valid.frequency <- rep(TRUE, length(infos))
		major.frequency <- rep(as.double(NA), length(infos))
		logger.warning(c("Could not detect MAF data; allele frequency is not considered"))
	} else {
		major.frequency <- rep(as.double(NA), length(infos))
		frequencies <- strsplit(gsub(regex.frequency, "\\1", infos[i]), ",", fixed = TRUE)
		major.frequency[i] <- sapply(frequencies, function(x) { max(as.double(x)) })
		rm(frequencies)
		is.valid.frequency <- (major.frequency <= 0.95)
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
	result$chromosome <- factor(result$chromosome, levels = CHROMOSOMES[[assembly]])
	result <- result[with(result, order(chromosome, location)), ]
	attr(result, "version") <- version.string
	result
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
		snps <- lapply(fnames, rnb.update.load.vcf, assembly = assembly)
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
	chromosomes <- names(snps)
	snps <- foreach(snp.df = snps) %dopar% rnb.construct.snp.types(snp.df)
	names(snps) <- chromosomes
	attr(snps, "version") <- db.version
	logger.status("Constructed tables of polymorphism types")

	return(snps)
}


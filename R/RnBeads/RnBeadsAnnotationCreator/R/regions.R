########################################################################################################################
## regions.R
## created: 2012-10-25
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Annotation of tiling regions, genes, promoters and other regulatory elements in the genome.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' download.ensembl.table
#'
#' Downloads a table from Biomart.
#'
#' @param database.name    Name of database that contains the table. This database should be among the ones listed by
#'                         \code{\link{listMarts()}}.
#' @param dataset.name     Name of the table to extract.
#' @param required.columns Columns of the table to download, in the form of a \code{character} vector. If specified, the
#'                         names of this vector are used to rename the columns of the downloaded table.
#' @param ...              Additional arguments passed to \code{\link{listMarts}} and \code{\link{useMart}}.
#' @return \code{data.frame} with the given columns, as returned by \code{\link{getBM}}.
#' @seealso \code{\link{useMart}} and \code{\link{getBM}} from package \pkg{biomaRt}
#'
#' @author Yassen Assenov
#' @noRd
download.ensembl.table <- function(database.name, dataset.name, required.columns, ...) {
	db.version <- listMarts(...)
	db.version <- as.character(db.version[db.version[, "biomart"] == database.name, "version"])
	mart <- useMart(database.name, dataset = dataset.name, ...)
	if (!all(required.columns %in% c(listAttributes(mart)[, "name"]))) {
		stop("Not all required attributes found in Ensembl")
	}
	result <- getBM(attributes = required.columns, mart = mart)
	if (!is.null(names(required.columns))) {
		colnames(result) <- names(required.columns)
	}
	attr(result, "version") <- db.version
	return(result)
}

########################################################################################################################

#' rnb.update.region.annotation.genes
#'
#' Creates a genomic annotation for gene and gene promoter regions.
#'
#' @param biomart.parameters Parameters passt to \code{\link{download.ensembl.table}} as a \code{list} containing at
#'                           least the following three items: \code{"database.name"}, \code{"dataset.name"} and
#'                           \code{"required.columns"}.
#' @return List of two items (\code{"genes"} and \code{"promoters"}), each of them is a list of containing the following
#'         two entries:
#'         \describe{
#'           \item{\code{"regions"}}{Object of type \code{GRangesList} containing \code{GRanges} objects for each
#'                chromosome. These ranges define region coordinates and annotation.}
#'           \item{\code{"mappings"}}{List of mappings for every element in \code{sites}. One mapping is a list of
#'                tables as objects of type \code{IRanges}, one per chromosome. Every table stores the range of indices 
#'                of \code{sites} on the respective chromosome that are contained in the corresponding region. Regions
#'                that do not contain sites are left out of the mappings.}.
#'         }
#'
#' @author Fabian Mueller
#' @noRd
rnb.update.region.annotation.genes <- function(biomart.parameters) {
	if (!suppressPackageStartupMessages(require(biomaRt))) {
		logger.error("Missing required package biomaRt")
	}
	get.single <- function(x, column.name) {
		if (length(unique(x)) != 1) {
			logger.error(c("Differing", column.name, "for Ensembl gene"))
		}
		return(x[1])
	}
	combine.multiple <- function(x) {
		if (all(is.na(x))) {
			return(NA)
		}
		return(paste(unique(sort(x[!is.na(x)])), collapse = ";"))
	}
	ensembl.genes <- do.call(download.ensembl.table, biomart.parameters)
	db.version <- attr(ensembl.genes, "version")
	ensembl.genes <- ensembl.genes[with(ensembl.genes, order(chromosome, start, end)), ]
	logger.status("Downloaded and sorted gene definitions from Ensembl using biomaRt")

	## Merge duplicate records
	ensembl.genes[ensembl.genes[, "symbol"] == "", "symbol"] <- NA
	gene.id.inds <- tapply(1:nrow(ensembl.genes), ensembl.genes[["id"]], identity)
	chroms <- sapply(gene.id.inds, function(i) { get.single(ensembl.genes[i, "chromosome"], "chromosomes") })
	chroms <- paste0("chr", chroms)
	starts <- sapply(gene.id.inds, function(i) { get.single(ensembl.genes[i, "start"], "start positions") })
	ends <- sapply(gene.id.inds, function(i) { get.single(ensembl.genes[i, "end"], "end positions") })
	strands <- sapply(gene.id.inds, function(i) { get.single(ensembl.genes[i, "strand"], "strands") })
	strands <- rnb.fix.strand(strands)
	symbols <- sapply(gene.id.inds, function(i) { combine.multiple(ensembl.genes[i, "symbol"]) })
	entrezids <- sapply(gene.id.inds, function(i) { combine.multiple(ensembl.genes[i, "entrezID"]) })
	CHROMOSOMES <- names(.globals[['CHROMOSOMES']])
	ensembl.genes.gr <- lapply(CHROMOSOMES, function(chrom) {
			cinds <- which(chroms == chrom)
			genes <- GRanges(seqnames = chrom,
				ranges = IRanges(start = starts[cinds], end = ends[cinds], names = names(gene.id.inds)[cinds]),
				strand = strands[cinds], symbol = symbols[cinds], entrezID = entrezids[cinds])
			seqlevels(genes) <- CHROMOSOMES
			seqlengths(genes) <- as.integer(seqlengths(rnb.genome.data())[CHROMOSOMES])
			genes <- rnb.sort.regions(genes)
			return(genes)
		}
	)
	names(ensembl.genes.gr) <- CHROMOSOMES
	rm(get.single, combine.multiple, ensembl.genes, gene.id.inds, chroms, starts, ends, strands, symbols, entrezids)
	ensembl.genes.gr <- GRangesList(ensembl.genes.gr)
	logger.status("Basic gene annotation completed")

	## Define promoters w.r.t. TSS
	## note the following property of TSS from the biomaRt ensembl annotations:
	## start_position is identical with min(transcript_start) and end_position is identical with max(transcript_end)
	## the TSS is transcript_start on the + strand and transcript_end on the - strand, respectively
	LENGTH.PROMOTER <- LENGTH.PROMOTER.UPSTREAM + LENGTH.PROMOTER.DOWNSTREAM
	ensembl.promoters.gr <- endoapply(ensembl.genes.gr, function(x) {
			resize(flank(x, LENGTH.PROMOTER.UPSTREAM), LENGTH.PROMOTER, fix = "start")
		}
	)
	logger.status("Basic promoter annotation completed")

	genome.data <- rnb.genome.data()
	ensembl.genes.gr <- append.cpg.stats(genome.data, ensembl.genes.gr)
	attr(ensembl.genes.gr, "version") <- db.version
	ensembl.promoters.gr <- append.cpg.stats(genome.data, ensembl.promoters.gr)
	attr(ensembl.promoters.gr, "version") <- db.version
	logger.status("Appended CpG-related statistics")

	return(list("genes" = ensembl.genes.gr, "promoters" = ensembl.promoters.gr))
}

########################################################################################################################

## rnb.update.region.annotation.tiling
##
## Creates a genomic annotation for tiling regions.
##
## @param assembly Genome assembly of interest.
## @return List with two entries:
##         \describe{
##           \item{\code{"regions"}}{Object of type \code{GRangesList} containing \code{GRanges} objects for each
##                chromosome. These ranges define region coordinates and annotation.}
##           \item{\code{"mappings"}}{List of mappings for every element in \code{sites}. One mapping is a list of
##                tables as objects of type \code{IRanges}, one per chromosome. Every table stores the range of indices 
##                of \code{sites} on the respective chromosome that are contained in the corresponding region. Regions
##                that do not contain sites are left out of the mappings.}.
##         }
##
## @author Fabian Mueller
rnb.update.region.annotation.tiling <- function(window.size=LENGTH.TILING){
	genome.data <- rnb.genome.data()
	CHROMOSOMES <- names(.globals[['CHROMOSOMES']])
	tiling.chrom <- function(chrom) {
		chromNames.gd <- match.chrom.names(CHROMOSOMES,seqnames(genome.data))
		chrom.length <- seqlengths(genome.data)[chromNames.gd[chrom]]
		seq.lengths <- seqlengths(genome.data)[chromNames.gd[CHROMOSOMES]]
		names(seq.lengths) <- CHROMOSOMES
		starts <- seq(1L, chrom.length, window.size)
		ends <- pmin(starts + window.size - 1L, chrom.length)
		GRanges(seqnames = chrom, ranges = IRanges(starts, ends),
			seqlengths = seq.lengths)
	}
	tiling.gr <- suppressWarnings(
		foreach(chrom = CHROMOSOMES, .packages = "GenomicRanges", .export = "CHROMOSOMES") %dopar%
			tiling.chrom(chrom))
	names(tiling.gr) <- CHROMOSOMES
	logger.status("Defined tiling regions for all supported chromosomes")
	tiling.gr <- append.cpg.stats(genome.data, tiling.gr)
	logger.status("Added CpG-related statistics to the tiling region definitions")
	return(tiling.gr)
}

########################################################################################################################

## rnb.update.download.cgis
##
## Downloads a CpG island definition table from the UCSC genome browser.
##
## @param download.url URL of the CGI table to download.
## @return CpG islands as an object of type \code{GRangesList}. Every element in this list corresponds to a chromosome.
##
## @author Yassen Assenov
rnb.update.download.cgis <- function(download.url) {
	genome.data <- rnb.genome.data()

	## Download the corresponding track from the UCSC Genome Browser
	cgis.file <- file.path(.globals[["DIR.PACKAGE"]], "temp", "cgis.txt.gz")
	if (file.exists(cgis.file)) {
		logger.info(c("File", cgis.file, "found; skipping download"))
	} else {
		if (download.file(download.url, cgis.file, quiet = TRUE, mode = "wb") != 0) {
			logger.error(c("Could not download", download.url))
		}
		logger.status(c("Downloaded CpG island track from", download.url))
	}

	## Post-process the table
	cgis.table <- read.delim(cgis.file, header = FALSE, stringsAsFactors = FALSE)
	cgis.table <- cgis.table[, c(2, 3, 4)] # focus only on chromosome, start and end columns
	cgis.table[, 2] <- cgis.table[, 2] + 1L # adjust start locations to 1-based

	## Convert the table to GRangesList
	CHROMOSOMES <- names(.globals[['CHROMOSOMES']])
	cgis.gr <- data.frame2GRanges(cgis.table, NULL, 1, 2, 3, NULL, assembly=NULL)
	cgis.gr <- cgis.gr[as.character(seqnames(cgis.gr)) %in% CHROMOSOMES]
	seqlevels(cgis.gr) <- CHROMOSOMES
	cgis.gr <- GenomicRanges::split(cgis.gr,seqnames(cgis.gr))
	chromNames.gd <- match.chrom.names(CHROMOSOMES,seqnames(genome.data))
	seqlengths(cgis.gr) <- as.integer(seqlengths(genome.data)[chromNames.gd[CHROMOSOMES]])
	cgis.gr <- append.cpg.stats(genome.data, rnb.sort.regions(cgis.gr))
	return(cgis.gr)
}

########################################################################################################################

#' rnb.update.region.annotation
#'
#' Initializes all region annotations for the given assembly.
#'
#' @param assembly Genome assembly to use.
#' @return List of genome region definitions. Every element in the list is dedicated to one region type and is an
#'         instance of \code{GRangesList}.
#'
#' @author Fabian Mueller
#' @noRd
rnb.update.region.annotation <- function(biomart.parameters, cgi.download.url =
		paste0(UCSC.FTP.BASE, .globals[['assembly']], "/database/cpgIslandExt.txt.gz")) {

	logger.start("Tiling Region Annotation")
	result <- list("tiling" = rnb.update.region.annotation.tiling())
	logger.completed()

	logger.start("Gene Annotation")
	result <- c(result, rnb.update.region.annotation.genes(biomart.parameters))
	logger.completed()

	logger.start("CpG Island Region Annotation")
	result[["cpgislands"]] <- rnb.update.download.cgis(cgi.download.url)
	logger.completed()

	attr(result, "builtin") <- sapply(result, function(x) { TRUE })
	result <- rnb.add.descriptions(result)
	return(result)
}

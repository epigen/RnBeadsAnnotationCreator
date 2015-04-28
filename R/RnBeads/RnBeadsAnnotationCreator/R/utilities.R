########################################################################################################################
## utilities.R
## created: 2014-02-13
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Helper functions used in constructing RnBeads annotation packages.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.genome.data
#'
#' Gets the targeted genome assembly sequence.
#'
#' @return Sequence data object for the current assembly.
#'
#' @author Yassen Assenov
#' @noRd
rnb.genome.data <- function() {
	if (is.null(.globals[['GENOME']])) {
		stop("undefined assembly")
	}
	get(.globals[['GENOME']])
}

########################################################################################################################

#' rnb.fix.strand
#'
#' Converts, if necessary, the provided strand information into a factor vector with levels \code{"+"}, \code{"-"} and
#' \code{"*"}.
#'
#' @param values Vector of strand information to be converted.
#' @return The (possibly modified) strand information as a factor vector.
#'
#' @author Yassen Assenov
#' @noRd
rnb.fix.strand <- function(values) {
	if (is.factor(values) && setequal(levels(values), c("+", "-", "*"))) {
		values[is.na(values)] <- "*"
	} else {
		values <- as.character(values)
		values[is.na(values)] <- "*"
		i.positive <- values %in% c("+", "1", "+1", "F", "f", "TOP")
		i.negative <- values %in% c("-", "-1", "R", "r", "BOT")
		values[i.positive] <- "+"
		values[i.negative] <- "-"
		values[!(i.positive | i.negative)] <- "*"
		values <- factor(values, levels = c("+", "-", "*"))
	}
	return(values)
}

########################################################################################################################

#' rnb.sort.regions
#'
#' Sorts the given regions based on start and end position.
#'
#' @param x Genomic regions as an object of type \code{\link{GRanges}} or \code{\link{GRangesList}}.
#' @return Set of the same regions as x, sorted based on start and end positions.
#'
#' @author Yassen Assenov
#' @noRd
rnb.sort.regions <- function(x) {
	if (inherits(x, "GRanges")) {
		return(x[order(start(x), end(x), as.integer(strand(x))), ])
	}
	if (inherits(x, "GRangesList")) {
		return(endoapply(x, function(y) { y[order(start(y), end(y), as.integer(strand(y))), ] }))
	}
	stop("invalid value for x")
}

########################################################################################################################

#' rnb.load.bed
#'
#' Loads a BED file into a \code{data.frame} with fixed column names. The file contents is validated for structure
#' (at least 3 columns), as well as for integer values in columns 2 and 3.
#'
#' @param fname BED file to load.
#' @return \code{data.frame} with at least 3 and at most 6 columns. The column names are: \code{"chromosome"},
#'         \code{"start"}, \code{"end"}, \code{"id"}, \code{"score"} and \code{"strand"}. Columns after the sixth one,
#'         if present, are dropped.
#'
#' @author Yassen Assenov
#' @noRd
rnb.load.bed <- function(fname) {
	BED.COLUMNS <- c("chromosome", "start", "end", "id", "score", "strand")
	tbl <- tryCatch(suppressWarnings(read.delim(fname, header = FALSE, quote = "", comment.char = "#",
				stringsAsFactors = FALSE, na.strings = "")), error = function(e) { e })
	if (inherits(tbl, "error")) {
		if (grepl("cannot open", tbl, fixed = TRUE)) {
			stop("cannot open file")
		}
		stop("invalid file format")
	}
	if (ncol(tbl) < 3) {
		stop("invalid file format; expected at least 3 columns")
	}
	if (!(is.integer(tbl[[2]]) && is.integer(tbl[[3]]))) {
		stop("invalid file format; expected start and end positions in columns 2 and 3, respectively")
	}
	if (length(BED.COLUMNS) < ncol(tbl)) {
		tbl <- tbl[, 1:length(BED.COLUMNS)]
	}
	colnames(tbl) <- BED.COLUMNS[1:ncol(tbl)]
	tbl[[1]] <- as.character(tbl[[1]])
	if (4 <= ncol(tbl)) {
		tbl[[4]] <- as.character(tbl[[4]])
		if (anyDuplicated(tbl[[4]]) == 0) {
			rownames(tbl) <- tbl[[4]]
			tbl <- tbl[, -4]
		}
	}
	return(tbl)
}

########################################################################################################################

#' get.cpg.stats
#'
#' Computes CpG-related statistics for the specified regions.
#'
#' @param chrom.sequence Chromosome sequence, usually obtained from the assembly's genome definition. This must be an
#'                       object of type \code{MaskedDNAString}.
#' @param starts         \code{integer} vector of start positions for the regions of interest. 
#' @param ends           \code{integer} vector of end positions for the regions of interest.
#' @return Table of statistics for the regions in the form of a \code{matrix} with the following columns:
#'         \code{"CpG"} and \code{"GC"}. The columns contain the number of CpG dinucleoties and the number of C and G
#'         bases in each region.
#'
#' @author Yassen Assenov
#' @noRd
get.cpg.stats <- function(chrom.sequence, starts, ends) {
	if (!(inherits(chrom.sequence, "MaskedDNAString") || inherits(chrom.sequence, "DNAString"))) {
		stop("invalid value for chrom.sequence")
	}
	if (!(is.integer(starts) && !any(is.na(starts)))) {
		stop("invalid value for starts")
	}
	if (!(is.integer(ends) && !any(is.na(ends)))) {
		stop("invalid value for ends")
	}
	chrom.regions <- suppressWarnings(Views(chrom.sequence, start = starts, end = ends))
	cg.freqs <- letterFrequency(chrom.regions, c("C", "G"))
	cbind(
		"CpG" = dinucleotideFrequency(chrom.regions)[, "CG"],
		"GC" = as.integer(rowSums(cg.freqs)),
		"C" = as.integer(cg.freqs[,"C"]),
		"G" = as.integer(cg.freqs[,"G"])
	)
}

########################################################################################################################

#' append.cpg.stats
#'
#' Appends additional metadata columns for CpG count and GC density to the specified regions.
#'
#' @param genome.data Genome of interest.
#' @param regionlist  Genomic regions as a list of \code{GRanges} objects (or an obect of type \code{GRangesList}),
#'                    containing one set of regions per chromosome.
#' @return The modified \code{regionlist}. Two columns are appedned to the metadata of each element in this list -
#'         \code{"CpG"} and \code{"GC"}. If the metadata already contains these columns, this function appends columns
#'         with similar names.
#'
#' @author Yassen Assenov
#' @noRd
append.cpg.stats <- function(genome.data, regionlist) {
	cpg.stats <- function(chrom) {
		chromName.gd <- match.chrom.names(chrom,seqnames(genome.data))
		stats <- get.cpg.stats(genome.data[[chromName.gd]], start(regionlist[[chrom]]), end(regionlist[[chrom]]))
		result <- regionlist[[chrom]]
		mcols(result) <- IRanges::cbind(mcols(result), DataFrame(stats))
		result
	}
	regions.enriched <- suppressWarnings(
		foreach(chrom = names(regionlist), .packages = "GenomicRanges",
			.export = c("get.cpg.stats", "genome.data", "regionlist")) %dopar% cpg.stats(chrom))
	names(regions.enriched) <- names(regionlist)
	return(GRangesList(regions.enriched))
}

########################################################################################################################

#' match.chrom.names
#'
#' invariant chromosome name matching. I.e. find the occurrences of the names in the first
#' vector in the second vector, neglecting the "chr" prefix
#'
#' @param q chromosome name query
#' @param q chromosome name reference
#' @return the names from r corresponding to the entries in q
#'
#' @author Fabian Mueller
#' @noRd
match.chrom.names <- function(q, r) {
	rr <- gsub("^chr","",r)
	qq <- gsub("^chr","",q)
	res <- r[match(qq,rr)]
	names(res) <- q
	return(res)
}

########################################################################################################################

#' createPackageScaffold
#' 
#' Creates a scaffold folder structure for an R package.
#'
#' @param pkg.name Name of the package to be created.
#' @param desc     Content of the DESCRIPTION file. This must be a named character vector with the headers as names.
#' @param dest     Destination directory where the package should be created.
#' @return Invisibly, \code{TRUE} if the package directory and its \code{DESCRIPTION} file were successfully created;
#'         \code{FALSE} otherwise.
#'
#' @author Fabian Mueller
#' @examples
#' createPackageScaffold("myPkg")
#' @noRd
createPackageScaffold <- function(
		pkg.name,
		desc=c(
			Package=pkg.name,
			Title=pkg.name,
			Description="automatically generated package",
			Author="Package Creator",
			Date=format(Sys.Date(), format="%Y-%m-%d"),
			License="GPL-3",
			Version="0.1"
		),
		dest=getwd()){
	pkg.base.dir <- file.path(dest,pkg.name)
	if (file.exists(pkg.base.dir)){
		stop("Package directory already exists")
	}
	## Create the folder structure
	for (dname in c("R", "data", "inst", "man", "temp")) {
		if (!dir.create(file.path(pkg.base.dir, dname), showWarnings = FALSE, recursive = TRUE)) {
			return(invisible(FALSE))
		}
	}
	## Create DESCRIPTION file
	desc.lines <- paste(names(desc),desc,sep=": ")
	writeLines(desc.lines,file.path(pkg.base.dir,"DESCRIPTION"))
	## Create NAMESPACE file
	writeLines(c(""),file.path(pkg.base.dir,"NAMESPACE"))
	## Create documentation (Rd) file
	assembly <- .globals[["assembly"]]
	fname <- system.file("extdata/templateRd.txt", package = "RnBeadsAnnotationCreator")
	txt <- gsub("%s", assembly, scan(fname, "", sep = "\n", na.strings = character(), quiet = TRUE), fixed = TRUE)
	fname <- paste0(pkg.base.dir, "/man/", assembly, ".Rd")
	cat(txt, file = fname, sep = "\n")
	## Create NEWS file
	fname <- system.file(paste0("extdata/NEWS.", assembly), package = "RnBeadsAnnotationCreator")
	if (file.exists(fname)) {
		txt <- scan(fname, "", sep = "\n", na.strings = character(), quiet = TRUE, blank.lines.skip = FALSE)
	} else {
		txt <- paste0("RnBeads.", assembly, " ", desc["Version"])
		txt <- c(txt, paste(rep("=", nchar(txt)), collapse = ""))
		txt <- c(txt, "",  paste0("* Initial release of RnBeads.", assembly, "."))
	}
	cat(txt, file = paste0(pkg.base.dir, "/inst/NEWS"), sep = "\n")
	invisible(TRUE)
}

########################################################################################################################

#' rnb.add.descriptions
#'
#' Adds an attribute \code{"description"} to each of the given annotation, if applicable.
#'
#' @param anns List of site or region annotation objects (instances of \code{GRangesList}).
#' @return The updated \code{list} of annotations with added or replaced attribute named \code{"description"}.
#'
#' @author Yassen Assenov
#' @noRd
rnb.add.descriptions <- function(anns) {
	for (aname in intersect(names(anns), names(ANNOT.DESCRIPTIONS))) {
		txt <- unname(ANNOT.DESCRIPTIONS[aname])
		txt.version <- attr(anns[[aname]], "version")
		if (isTRUE(txt.version != "")) {
			txt <- paste0(txt, ", version ", txt.version)
		}
		attr(anns[[aname]], "description") <- txt
	}
	anns
}

########################################################################################################################

#' rnb.create.mappings
#' 
#' Creates mappings from every possible region annotation to every possible site annotation.
#'
#' @param regions \code{list} of region annotations; every annotation should be a \code{GRangesList} object storing one
#'                \code{GRanges} instance per chromosome.
#' @param sites   \code{list} of site annotations; every annotation should be a \code{GRangesList} object storing one
#'                \code{GRanges} instance per chromosome.
#' @return The initialized mapping structure, invisibly.
#'
#' @author Yassen Assenov
#' @noRd
rnb.create.mappings <- function(regions = .globals[['regions']], sites = .globals[['sites']]) {
	mappings <- 
		suppressWarnings(foreach(re = regions) %:% foreach(si = sites) %dopar% RnBeads:::rnb.regions2sites(re, si))
	names(mappings) <- names(regions)
	for (rname in names(mappings)) {
		names(mappings[[rname]]) <- names(sites)
	}
	assign('mappings', mappings, .globals)
	return(invisible(mappings))
}

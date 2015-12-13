########################################################################################################################
## createAnnotationPackage.R
## created: 2014-02-13
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Creation of annotation packages for RnBeads.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' update.annot
#'
#' Updates an annotation table.
#'
#' @param object.name Name of the annotation object to load (if it exists) or save (if it needs to be created).
#' @param info        Human-readable description of the annotation. This must be a one-element \code{character} vector.
#' @param unpdate.fun Function to be called for creating the annotation object.
#' @param ...         Parameters passed to the updating function.
#' @return The loaded or initialized annotation object, invisibly.
#'
#' @author Yassen Assenov
#' @noRd
update.annot <- function(object.name, info, update.fun, ...) {
	location <- file.path(.globals[['DIR.PACKAGE']], "temp")
	fname <- file.path(location, paste0(object.name, ".RDS"))
	if (file.exists(fname)) {
		obj <- readRDS(fname)
		logger.status(c("Loaded", info, "from", fname))
	} else {
		obj <- update.fun(...)
		con <- gzfile(fname, "wb", compression = 9L)
		saveRDS(obj, file = con)
		close(con)
		logger.status(c("Saved", info, "to", fname))
	}
	assign(object.name, obj, .globals)
	return(invisible(obj))
}

########################################################################################################################

#' rnb.get.package.data.file
#'
#' Gets the full name of a data file (dedicated to a specific site or region annotation) in the currently processed
#' annotation package.
#'
#' @param annotation.name Name of the major annotation table.
#' @return Full path to the \code{.RData} file that stores or should store the specified annotation table.
#'
#' @author Yassen Assenov
#' @noRd
rnb.get.package.data.file <- function(annotation.name) {
	file.path(.globals[['DIR.PACKAGE']], "data", paste0(annotation.name, ".RData"))
}

########################################################################################################################

#' rnb.export.annotations.to.data.files
#'
#' Exports the annotation tables and mappings in the environment \code{.globals} to the \code{"data"} directory of the
#' annotation package currently being constructed.
#'
#' @author Yassen Assenov
#' @noRd
rnb.export.annotations.to.data.files <- function() {
	sites.full <- .globals[['sites']]
	framework <- list(
		"GENOME" = .globals[['GENOME']],
		"CHROMOSOMES" = .globals[['CHROMOSOMES']],
		"regions" = list(),
		"sites" = lapply(sites.full, function(x) { NULL }),
		"controls" = list(),
		"mappings" = list())
	for (sname in names(sites.full)) {
		sites <- list("sites" = sites.full[[sname]], "mappings" = lapply(.globals[['mappings']], "[[", sname))
		if (grepl("^probes", sname)) {
			platform.id <- gsub("^probes", "", sname)
			sites[[paste0("controls", platform.id)]] <- get(sname, .globals)[["controls"]]
			framework[["controls"]][paste0("controls", platform.id)] <- list(NULL)
		}
		save(sites, file = rnb.get.package.data.file(sname), compression_level = 9L)
		logger.status(c("Saved", sname, "annotation table and mappings to the package data"))
	}
	if (length(framework[["controls"]]) == 0) {
		framework[["controls"]] <- NULL
	} else {
		cinfo <- sub("^controls", "probes", names(framework[["controls"]]))
		names(cinfo) <- names(framework[["controls"]])
		attr(framework[["controls"]], "sites") <- cinfo
		rm(cinfo)
	}
	assembly <- .globals[['assembly']]
	.globals[[assembly]] <- framework
	save(list = assembly, file = rnb.get.package.data.file(assembly), envir = .globals, compression_level = 9L)
	regions <- .globals[['regions']]
	save(regions, file = rnb.get.package.data.file("regions"), compression_level = 9L)
	logger.status("Saved region annotation table to the package data")
}

########################################################################################################################

#' createAnnotationPackage
#'
#' Generates an R package containing the RnBeads annotation for the specified genome assembly.
#'
#' @param assembly Targeted genome assembly. Must be one of \code{"hg38"}, \code{"hg19"}, \code{"mm10"}, \code{"mm9"},
#'                 \code{"rn5"}.
#' @param dest     Destination directory where the package should be generated.
#' @param maxSNPs  Maximum number of dbSNP records in preprocessed table for a single chromosome, given as a single
#'                 non-negative \code{integer} value. If a larger table (than this threshold) is found, the annotation
#'                 creator stops and suggests a solution that utilizes the environment of a computational cluster.
#'                 Setting this parameter to \code{0} disables such a test. This parameter is used only when SNP data is
#'                 available for the targeted genome assembly.
#' @param cleanUp  Flag indicating if the temporary directory in the generated package (containing downloaded and
#'                 partially processed resource files) is to be removed after all annotation structures are created.
#' @return None (invisible \code{NULL}).
#' @author Fabian Mueller
#' @examples
#' \donttest{
#' createAnnotationPackage("hg38")
#' }
#' @export
createAnnotationPackage <- function(assembly,dest=getwd(),maxSNPs=500000L,cleanUp=TRUE){
	## Validate parameters
	if (!(is.character(assembly) && length(assembly) == 1 && (isTRUE(assembly != "")))) {
		stop("invalid value for assembly")
	}
	if (!(is.character(dest) && length(dest) == 1 && (isTRUE(dest != "")))) {
		stop("invalid value for dest")
	}
	if (is.double(maxSNPs) && isTRUE(all(maxSNPs == as.integer(maxSNPs)))) {
		maxSNPs <- as.integer(maxSNPs)
	}
	if (!(is.integer(maxSNPs) && length(maxSNPs) == 1 && isTRUE(0 <= maxSNPs))) {
		stop("invalid value for maxSNPs")
	}
	if (!(is.logical(cleanUp) && length(cleanUp) == 1 && (!is.na(cleanUp)))) {
		stop("invalid value for cleanUp")
	}

	## Validate the genome assembly
	function.name <- paste0("createAnnotationPackage.", assembly)
	if (!exists(function.name)) {
		stop("Unsupported assembly")
	}

	## Initialize the environment variables
	rm(list = ls(.globals), pos = .globals)
	assign('assembly', assembly, .globals)
	dir.package <- file.path(dest, paste0("RnBeads.", assembly))
	assign('DIR.PACKAGE', dir.package, .globals)
	assign('SNP_MAX', maxSNPs, .globals)

	## Initialize package directories
	if (file.exists(dir.package)) {
		dir.package.state <- "(existing)"
	} else {
		pkgName <- paste0("RnBeads.", assembly)
		txt <- paste0("Automatically generated RnBeads annotation package for the assembly ", assembly, ".")
		createdScaffold <- createPackageScaffold(
			pkgName,
			desc=c(
				Package=pkgName,
				Title=pkgName,
				Description=txt,
				Author="RnBeadsAnnotationCreator",
				Maintainer="RnBeadsAnnotationCreator <rnbeads@mpi-inf.mpg.de>",
				Date=format(Sys.Date(), format="%Y-%m-%d"),
				License="GPL-3",
				Version="0.1",
				Depends="\n\tR (>= 3.0.0),\tGenomicRanges",
				Suggests="\n\tRnBeads"
			),
			dest = dest
		)
		if (!createdScaffold) {
			stop("Could not create package directory")
		}
		dir.package.state <- "(created)"
	}

	## Create the annotation package
	logging2console <- (!logger.isinitialized())
	if (logging2console) {
		logger.start(fname = NA)
	}
	logger.start("Creating Annotation Package")
	logger.info(c("Assembly:", assembly))
	logger.info(c("Package directory", dir.package.state, ":", dir.package))
	rm(dir.package, dir.package.state)
	do.call(function.name, list())

	## Clean the temporary directory
	if (cleanUp){
		if (unlink(file.path(.globals[['DIR.PACKAGE']], "temp"), recursive = TRUE) != 0L) {
			logger.warning("Could not clean package temporary directory")
		} else {
			logger.status("Cleaned package temporary directory")
		}
	}
	logger.completed()
	if (logging2console) {
		logger.close()
	}
}

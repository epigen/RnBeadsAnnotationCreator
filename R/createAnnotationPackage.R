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
#' @return The loaded or initialized annotation object.
#' 
#' @author Yassen Assenov
#' @noRd
update.annot <- function(object.name, info, update.fun, ...) {
	location <- file.path(.globals[['DIR.PACKAGE']], 'temp')
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
	return(obj)
}

########################################################################################################################

#' createAnnotationPackage
#' 
#' Generates an R package containing the RnBeads annotation for the specified genome assembly.
#'
#' @param assembly    Targeted genome assembly. Must be one of \code{"hg38"}, \code{"hg19"}, \code{"mm10"},
#'                    \code{"mm9"}, \code{"rn5"}.
#' @param dest        Destination directory where the package should be generated.
#' @param cores.count Number of processing cores to be used in the computations.
#' @author Fabian Mueller
#' @examples
#' createAnnotationPackage("hg38")
#' @export
createAnnotationPackage <- function(assembly,dest=getwd(),cores.count=1L){
	## Validate parameters
	if (!(is.character(assembly) && length(assembly) == 1 && (isTRUE(assembly != "")))) {
		stop("invalid value for assembly")
	}
	if (!(is.character(dest) && length(dest) == 1 && (isTRUE(dest != "")))) {
		stop("invalid value for dest")
	}
	if (is.double(cores.count) && isTRUE(all(cores.count == as.integer(cores.count)))) {
		cores.count <- as.integer(cores.count)
	}
	if (!(is.integer(cores.count) && length(cores.count) == 1 && (isTRUE(cores.count >= 1L)))) {
		stop("invalid value for cores.count")
	}

	## Validate the genome assembly
	function.name <- paste0("createAnnotationPackage.", assembly)
	if (!exists(function.name)) {
		stop("Unsupported assembly")
	}

	## Initialize the environment variables and the logger
	assign("assembly", assembly, .globals)
	dir.package <- file.path(dest, paste0("RnBeads.", assembly))
	assign("DIR.PACKAGE", dir.package, .globals)

	## Initialize parallel processing
	if (cores.count != 1) {
		registerDoParallel(cores.count)
	}

	## Initialize package directory
	if (file.exists(dir.package)) {
		dir.package.state <- "(existing)"
	} else {
		if (!createPackageScaffold(paste0("RnBeads.", assembly), dest = dest)) {
			rm(list = c("assembly", "DIR.PACKAGE"), pos = .globals)
			stop("Could not create package directory")
		}
		dir.package.state <- "(created)"
	}

	## Create the annotation package
	logger.start("Creating Annotation Package", fname = file.path(dir.package, "temp", "annotation.log"))
	logger.info(c("Assembly:", assembly))
	logger.info(c("Package directory", dir.package.state, ":", dir.package))
	rm(dir.package, dir.package.state)
	do.call(function.name, list(dest = dest))
}

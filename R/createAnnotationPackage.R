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
#' @param object.name ...
#' @param info        Human-readable description of the annotation stored in a one-element \code{character} vector.
#' @param unpdate.fun Function to be called for creating the annotation object.
#' @param ...         Parameters to pass to the updating function.
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
#' @return invisible \code{TRUE} if successful
#' @author Fabian Mueller
#' @examples
#' createAnnotationPackage("hg38")
#' @export
createAnnotationPackage <- function(assembly,dest=getwd(),cores.count=1L){
	assign('assembly', assembly, .globals)
	assign('DIR.PACKAGE', file.path(paste0("RnBeads.", assembly), dest), .globals)
	logger.start('Creating Annotation Package', fname = NA)
	logger.info(c('Assembly:', assembly))
	if (cores.count != 1) {
		registerDoParallel(cores.count)
	}
	createPackageScaffold(paste0("RnBeads.", assembly), dest = dest)
	if (assembly == "hg38") {
		createAnnotationPackage.hg38(dest)
	}
	logger.completed()
	invisible(TRUE)
}

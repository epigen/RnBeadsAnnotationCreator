#' createAnnotationPackage
#' 
#' Given an assembly identifier, this function generates an R-package containing the
#' RnBeads annotation for that assembly
#' @param assembly identifier for the genome assembly
#' @param dest destination directory where the package should be generated
#' @return invisible \code{TRUE} if successful
#' @author Yassen Assenov, Fabian Mueller
#' @export 
#' @examples 
#' createAnnotationPackage("hg38")
createAnnotationPackage <- function(assembly,dest=getwd()){
	if (assembly == "hg38"){
		createAnnotationPackage.hg38(dest)
	}
	invisible(TRUE)
}

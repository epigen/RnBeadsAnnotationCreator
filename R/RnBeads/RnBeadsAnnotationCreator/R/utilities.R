#' createPackageScaffold
#' 
#' Creates a scaffold folder structure for an R-package
#' @param pkg.name Name of the package to be created
#' @param desc Content of the DESCRIPTION file. Should be a named character vector with the headers as names.
#' @param dest destination directory where the package should be created
#' @return invisible \code{TRUE} if successful
#' @author Fabian Mueller
#' @examples 
#' createPackageScaffold("myPkg")
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
		stop("Package already exists")
	}
	#create folder structure
	dir.create(pkg.base.dir)
	dir.create(file.path(pkg.base.dir,"R"))
	dir.create(file.path(pkg.base.dir,"man"))
	dir.create(file.path(pkg.base.dir,"inst"))
	dir.create(file.path(pkg.base.dir,"data"))
	#create the DESCRIPTION file
	desc.lines <- paste(names(desc),desc,sep=": ")
	writeLines(desc.lines,file.path(pkg.base.dir,"DESCRIPTION"))
	invisible(TRUE)
}

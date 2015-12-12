#' createAnnotationPackage
#'
#' Generates an R package containing the RnBeads annotation for the specified genome assembly.
#'
#' @param example Name of the example dataset to be created. currently supported datasets: "ziller2011_450K"
#' @param dest     Destination directory where the package should be generated.
#' @return the example RnBSet object.
#' @author Fabian Mueller
#' @examples
#' \donttest{
#' createExampleDataset()
#' }
#' @export
createExampleDataset <- function(example="ziller2011_450K", dest=getwd()){
	rnbs <- NULL
	if (example=="ziller2011_450K"){
		dataset.url <- "http://rnbeads.mpi-inf.mpg.de/publication/data/ziller21011_450K_subset.RData"
		#TODO: read dataset from url
		downloaded.fn <- tempfile(pattern="exampleDataset", fileext=".RData")
		download.file(dataset.url, destfile=downloaded.fn)
		load(downloaded.fn)
		rnbs.orig <- rnb.set.example

		aa <- annotation(rnbs.orig)
		rnbs <- RnBeadRawSet(
			pheno(rnbs.orig),
			probes=rownames(aa),
			M=rnbs.orig@M,
			U=rnbs.orig@U,
			M0=rnbs.orig@M0,
			U0=rnbs.orig@U0,
			bead.counts.M=rnbs.orig@bead.counts.M,
			bead.counts.U=rnbs.orig@bead.counts.U,
			p.values=rnbs.orig@pval.sites,
			qc = rnbs.orig@qc,
			platform = "450k",
			beta.offset=100,
			summarize.bead.counts=TRUE,
			summarize.regions=TRUE,
			region.types = rnb.region.types.for.analysis("hg19"),
			useff=FALSE,
			ffcleanup=FALSE
		)

	} else {
		stop("Unknown example dataset identifier")
	}
	rnb.set.example <- rnbs
	if (!is.null(rnb.set.example)) {
		save(rnb.set.example, file=file.path(dest,paste0(example,".RData")), compression_level = 9L)
	}
	return(rnb.set.example)
}

########################################################################################################################
## sites.R
## created: 2012-02-13
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions dedicated to annotating loci relevant for bisulfite sequencing.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

## rnb.update.sites
##
## Creates lists of Genomic ranges for genomic sites targeted by methylation.
##
## @param cpgislands Region annotation of the CpG islands. If this is specified, the sites annotation is enriched with
##                   a column named \code{"CGI Relation"}.
## @return List of \code{GRangesList} objects. Every item is dedicated to a motif (e.g. CpGs), and groups \code{GRanges}
##         instances, one per chromosome.
## @author Fabian Mueller
rnb.update.sites <- function(cpgislands = NULL) {
	genome.data <- get.genome.data()
	CHROMOSOMES <- .globals[["CHROMOSOMES"]]
	chrom.lengths <- seqlengths(genome.data)[CHROMOSOMES]
	sites <- list()
	pp.dnas <- DNAStringSet(NUCLEOTIDE.PATTERNS)
	for (i in names(NUCLEOTIDE.PATTERNS)){
		pp.p <- pp.dnas[[i]]
		pp.m <- reverseComplement(pp.p)
		curSites <- lapply(CHROMOSOMES, function(chrom) {
			matches.st <- lapply(list(pp.p, pp.m), function(x) {
					ranges(matchPattern(x, genome.data[[chrom]]))
			})
			cp.starts <- start(matches.st[[1]]) - LENGTH.NEIGHBORHOOD %/% 2L + 1L
			cp.ends <- cp.starts + LENGTH.NEIGHBORHOOD - 1L
			cpg.stats <- suppressWarnings(get.cpg.stats(genome.data[[chrom]], cp.starts, cp.ends))
			matches.gr <- mapply(GRanges, seqnames = list(chrom), ranges = matches.st, strand = list("+", "-"),
				"CpG" = list(cpg.stats[, "CpG"]), "GC" = list(cpg.stats[, "GC"]))
			matches.gr <- rnb.sort.regions(do.call(c, matches.gr))
#			names(matches.gr) <- paste("cg", chrom, start(matches.gr), strand(matches.gr), sep = ".")
			seqlevels(matches.gr) <- names(CHROMOSOMES)
			seqlengths(matches.gr) <- chrom.lengths
			matches.gr
		})
		sites[[i]] <- GRangesList(curSites)
		logger.status(c("Created site annotation for", i))
	}

	if (!is.null(cpgislands)) {
		sites <- rnb.update.site.annotation.with.cgistatus(sites, cpgislands)
		logger.status("Enriched sites with CpG island information")
	}
	return(sites)
}

########################################################################################################################

## rnb.update.site.annotation.with.probes
##
## Enriches the annotation of genomic sites by adding columns showing if they are covered by the included platforms.
##
## @param sites          List of \code{GRangesList} objects, one per site type or probe annotation.
## @param query.probes   \code{character} vector storing the names of all elements in \code{sites} that define probe
##                       annotations.
## @param platform.names Platform names for the probe annotations. This must be a \code{character} vector of the same
##                       length as \code{query.probes}. This parameter effectively specifies the metadata columns in the
##                       sites annotation that will be added or replaced.
## @return The modified \code{sites}.
##
## @author Yassen Assenov
rnb.update.site.annotation.with.probes <- function(sites, query.probes = "probes450",
	platform.names = "HumanMethylation450") {

	for (i in 1:length(query.probes)) {
		qprobes <- query.probes[i]
		for (context.name in setdiff(names(sites), query.probes)) {
			chroms <- intersect(names(sites[[context.name]]), names(sites[[qprobes]]))
			overlaps <- mapply(GenomicRanges::findOverlaps, sites[[context.name]][chroms], sites[[qprobes]][chroms],
				MoreArgs = list(minoverlap = 2L), SIMPLIFY = FALSE)
			result <- list()
			for (chrom in names(sites[[context.name]])) {
				result[[chrom]] <- sites[[context.name]][[chrom]]
				platform.indices <- rep.int(as.integer(NA), length(result[[chrom]]))
				if (chrom %in% names(overlaps)) {
					platform.indices[queryHits(overlaps[[chrom]])] <- subjectHits(overlaps[[chrom]])
				}
				mcols(result[[chrom]])[, platform.names[i]] <- platform.indices
			}
			sites[[context.name]] <- GRangesList(result)
			rm(chroms, overlaps, result, chrom, platform.indices)
		}
	}
	return(sites)
}

########################################################################################################################

## rnb.update.site.annotation.with.cgistatus
##
## Enriches the annotation of genomic sites by adding columns showing their relation to CpG island status. This
## relation is one of \code{"Island"}, \code{"Shore"}, \code{"Shelf"}, \code{"Open Sea"}.
##
## @param sites      List of \code{GRangesList} objects, one per site type or probe annotation.
## @param cpgislands \code{GRangesList} object specifying the locations of CpG islands.
## @return The modified \code{sites}.
##
## @author Fabian Mueller
rnb.update.site.annotation.with.cgistatus <- function(sites, cpgislands) {
	sites.new <- sites
	chromosomes <- names(cpgislands)

	enrich.f <- function(ss, cgis) {
		cgi.relation <- rep.int("Open Sea", length(ss))
		if (length(cgis) != 0) {
			shores <- c(flank(cgis,LENGTH.CGI.SHORE-1L,start=TRUE),flank(cgis,LENGTH.CGI.SHORE-1L,start=FALSE))
			shelves <- c(flank(shores,LENGTH.CGI.SHELF-1L,start=TRUE),flank(shores,LENGTH.CGI.SHELF-1L,start=FALSE))
			#in the following order is essential
			is.shelf <- ss %in% shelves
			cgi.relation[is.shelf] <- "Shelf"
			is.shore <- ss %in% shores
			cgi.relation[is.shore] <- "Shore"
			is.cgi <- ss %in% cgis
			cgi.relation[is.cgi] <- "Island"
		}
		elementMetadata(ss)[, "CGI Relation"] <-
			factor(cgi.relation, levels = c("Open Sea", "Shelf", "Shore", "Island"))
		ss
	}

	for (i in names(sites)) {
		if (!setequal(chromosomes, names(sites[[i]]))) {
			stop("Incompatible sites and CGI ranges")
		}
		grl <- foreach(ss = as.list(sites[[i]][chromosomes]), cgis = as.list(cpgislands)) %dopar% enrich.f(ss, cgis)
		names(grl) <- chromosomes
		## Fix issue with column being renamed to CGI.Relation
		grl <- unlist(GRangesList(grl))
		i.cgirelation <- which(colnames(mcols(grl)) == "CGI.Relation")
		if (length(i.cgirelation) != 0) {
			colnames(mcols(grl))[i.cgirelation] <- "CGI Relation"
		}
		sites.new[[i]] <- GenomicRanges::split(grl, seqnames(grl))
	}
	
	return(sites.new)
}

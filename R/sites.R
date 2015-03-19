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
	chromNames.gd <- match.chrom.names(CHROMOSOMES,seqnames(genome.data))
	chrom.lengths <- seqlengths(genome.data)[chromNames.gd[CHROMOSOMES]]
	names(chrom.lengths) <- CHROMOSOMES
	sites <- list()
	pp.dnas <- DNAStringSet(NUCLEOTIDE.PATTERNS)
	for (i in names(NUCLEOTIDE.PATTERNS)){
		pp.p <- pp.dnas[[i]]
		pp.m <- reverseComplement(pp.p)
		curSites <- lapply(CHROMOSOMES, function(chrom) {
			matches.st <- lapply(list(pp.p, pp.m), function(x) {
					ranges(matchPattern(x, genome.data[[chromNames.gd[chrom]]]))
			})
			cp.starts <- start(matches.st[[1]]) - LENGTH.NEIGHBORHOOD %/% 2L + 1L
			cp.ends <- cp.starts + LENGTH.NEIGHBORHOOD - 1L
			cpg.stats <- suppressWarnings(get.cpg.stats(genome.data[[chromNames.gd[chrom]]], cp.starts, cp.ends))
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

	## TODO: Update sites with SNP information if present

	return(rnb.add.descriptions(sites))
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
	chromosomes <- names(cpgislands)

	enrich.f <- function(ss, cgis) {
		cgi.relations <- c("Open Sea", "Shelf", "Shore", "Island")
		cgi.relation <- factor(rep.int(cgi.relations[1], length(ss)), levels = cgi.relations)
		if (length(cgis) != 0) {
			shifts <- c(LENGTH.CGI.SHELF + LENGTH.CGI.SHORE, LENGTH.CGI.SHORE, 0L)
			names(shifts) <- cgi.relations[-1]
			x <- start(ss)
			for (j in names(shifts)) {
				subject <- GRanges(seqnames(cgis), IRanges(start(cgis) - shifts[j], end(cgis) - shifts[j]))
				cgi.relation[overlapsAny(ss, subject, ignore.strand = TRUE)] <- j
			}
		}
		elementMetadata(ss)[, "CGI Relation"] <- cgi.relation
		ss
	}

	sites.new <- list()
	for (i in names(sites)) {
		if (!setequal(chromosomes, names(sites[[i]]))) {
			stop("Incompatible sites and CGI ranges")
		}
		grl <- foreach(ss = as.list(sites[[i]][chromosomes]), cgis = as.list(cpgislands),
				.export = c("LENGTH.CGI.SHELF", "LENGTH.CGI.SHORE")) %dopar% enrich.f(ss, cgis)
		names(grl) <- chromosomes
		sites.new[[i]] <- GRangesList(grl)
	}
	
	return(sites.new)
}

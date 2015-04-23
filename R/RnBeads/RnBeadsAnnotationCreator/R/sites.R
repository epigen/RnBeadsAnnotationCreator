########################################################################################################################
## sites.R
## created: 2012-02-13
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions dedicated to annotating loci relevant for bisulfite sequencing.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.sites
#'
#' Creates lists of Genomic ranges for genomic sites targeted by methylation.
#'
#' @param cpgislands Region annotation of the CpG islands. If this is specified, the sites annotation is enriched with
#'                   a column named \code{"CGI Relation"}.
#' @param snps       SNP records as a \code{list} of \code{data.frame}s, one per chromosome.
#' @return List of \code{GRangesList} objects. Every item is dedicated to a motif (e.g. CpGs), and groups \code{GRanges}
#'         instances, one per chromosome.
#' @author Fabian Mueller
#' @noRd
rnb.update.sites <- function(cpgislands = .globals[['regions']][['cpgislands']], snps = .globals[['snps']]) {
	genome.data <- rnb.genome.data()
	CHROMOSOMES <- names(.globals[['CHROMOSOMES']])
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
			seqlevels(matches.gr) <- CHROMOSOMES
			seqlengths(matches.gr) <- chrom.lengths
			matches.gr
		})
		names(curSites) <- CHROMOSOMES
		sites[[i]] <- GRangesList(curSites)
		logger.status(c("Created site annotation for", i))
	}

	## Update sites with CGI annotation if present
	if (!is.null(cpgislands)) {
		sites <- rnb.update.site.annotation.with.cgistatus(sites, cpgislands)
		logger.status("Enriched sites with CpG island information")
	}

	## Update sites with SNP information if present
	if (!is.null(snps)) {
		sites <- rnb.update.site.annotation.with.snps(sites, snps)
	}

	rnb.add.descriptions(sites)
}

########################################################################################################################

#' rnb.update.site.annotation.with.probes
#'
#' Enriches the annotation of genomic sites by adding columns showing if they are covered by the included platforms.
#'
#' @param sites          List of \code{GRangesList} objects, one per site type or probe annotation.
#' @param query.probes   \code{character} vector storing the names of all elements in \code{sites} that define probe
#'                       annotations.
#' @param platform.names Platform names for the probe annotations. This must be a \code{character} vector of the same
#'                       length as \code{query.probes}. This parameter effectively specifies the metadata columns in the
#'                       sites annotation that will be added or replaced.
#' @return The modified \code{sites}.
#'
#' @author Yassen Assenov
#' @noRd
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
		logger.status(c("Computed overlaps with", qprobes))
	}
	return(sites)
}

########################################################################################################################

#' rnb.update.site.annotation.with.cgistatus
#'
#' Enriches the annotation of genomic sites by adding columns showing their relation to CpG island status. This
#' relation is one of \code{"Island"}, \code{"Shore"}, \code{"Shelf"}, \code{"Open Sea"}.
#'
#' @param sites      List of \code{GRangesList} objects, one per site type or probe annotation.
#' @param cpgislands \code{GRangesList} object specifying the locations of CpG islands.
#' @return The modified \code{sites}.
#'
#' @author Fabian Mueller
#' @noRd
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
		grl <- suppressWarnings(
			foreach(ss = as.list(sites[[i]][chromosomes]), cgis = as.list(cpgislands),
				.export = c("LENGTH.CGI.SHELF", "LENGTH.CGI.SHORE")) %dopar% enrich.f(ss, cgis))
		names(grl) <- chromosomes
		sites.new[[i]] <- GRangesList(grl)
	}
	
	return(sites.new)
}

########################################################################################################################

#' rnb.update.site.annotation.with.snps
#'
#' Enriches the annotation of genomic sites by adding a column showing if they overlap with records from dbSNP.
#'
#' @param sites List of \code{GRangesList} objects, one per site type or probe annotation.
#' @param snps  SNP records as a \code{list} of \code{data.frame}s, one per chromosome.
#' @return The modified \code{sites}.
#'
#' @details This function adds or replaces the annotation column \code{"SNP"} to every site annotation. For every CpG
#'          dinucleotide or probe, the SNP value is a \code{character} string containing a comma-separated lits of
#'          identifiers of dbSNP records (row names of the respective \code{data.frame} in \code{snps}) that overlap
#'          with the CpG.
#'
#' @author Yassen Assenov
#' @noRd
rnb.update.site.annotation.with.snps <- function(sites, snps) {
	chromosomes <- names(snps)

	for (i in names(sites)) {
		sites[[i]] <- endoapply(sites[[i]], function(x) { mcols(x)[, "SNPs"] <- as.character(NA); x })
		logger.status(c("Added column SNP to the annotation of sites", i))
		for (chrom in names(sites[[i]])) {
			cpg.coords <- start(sites[[i]][[chrom]])
			if (chrom %in% chromosomes) {
				cpg.coords <- cpg.coords[seq(1L, length(cpg.coords), by = 2L)]
				dfr <- unique(snps[[1]][, 1:2])

				hits <- findOverlaps(IRanges(cpg.coords, cpg.coords + 1L), IRanges(dfr[, 1], dfr[, 2]))
				lbls <- tapply(subjectHits(hits), queryHits(hits),
					function(j) { paste(rownames(dfr)[j], collapse = ",") })
				j <- rep(as.integer(names(lbls)), each = 2L) * 2L - c(1L, 0L)
				mcols(sites[[i]][[chrom]])[j, "SNPs"] <- rep(unname(lbls), each = 2L)
				rm(dfr, hits, lbls, j)
			}
		}
		logger.status(c("Overlapped", i, "sites with SNPs"))
		rm(chrom, cpg.coords)
	}
	return(sites)
} 

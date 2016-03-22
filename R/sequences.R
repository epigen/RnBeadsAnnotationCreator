########################################################################################################################
## sequences.R
## created: 2015-12-10
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Utility functions for working with genomic sequences encoded as character vectors.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

rnb.seq.bisulfite.convert <- function(x.char) {
	for (i in 1:length(x.char)) {
		xx <- x.char[[i]]
		for (j in 1:length(xx)) {
			if (xx[j] == "C" && (j == length(xx) || xx[j + 1] != "G")) {
				xx[j] <- "T"
			}
		}
		x.char[[i]] <- xx
	}
	x.char
}

########################################################################################################################

rnb.seq.reverse.complement <- function(x.char) {
	lapply(x.char, function(xx) { rev(unname(c("A" = "T", "C" = "G", "G" = "C", "T" = "A", "N" = "N")[xx])) })
}

########################################################################################################################

rnb.seq.replace <- function(x.char, s.old, s.new) {
	for (i in 1:length(x.char)) {
		xx <- x.char[[i]]
		xx[xx == s.old] <- s.new
		x.char[[i]] <- xx
	}
	x.char
}

########################################################################################################################

rnb.infinium.probe.coords <- function(loci, pr.design, pr.strand) {
	if (pr.design == "I") {
		if (pr.strand == "+") {
			return(cbind(loci + 0L, loci + 50L))
		}
		return(cbind(loci - 48L, loci + 2L))
	} # else pr.design == "II"
	if (pr.strand == "+") {
		return(cbind(loci + 1L, loci + 51L))
	}
	return(cbind(loci - 49L, loci + 1L))
}

########################################################################################################################

rnb.seq.get.expected.sequence <- function(chrom.sequence, loci, pr.design, pr.strand) {
	p.coords <- rnb.infinium.probe.coords(loci, pr.design, pr.strand)
	dna.seq <- suppressWarnings(Views(chrom.sequence, start = p.coords[, 1], end = p.coords[, 2]))
	dna.seq <- strsplit(as.character(dna.seq), NULL)
	if (pr.strand == "-") {
		dna.seq <- rnb.seq.reverse.complement(dna.seq)
	}
	alleles.exp <- rnb.seq.reverse.complement(rnb.seq.bisulfite.convert(dna.seq))
	if (pr.design == "I") {
		alleles.exp <- list("A" = rnb.seq.replace(alleles.exp, "G", "A"), "B" = alleles.exp)
	} else { # pr.design == "II"
		alleles.exp <- rnb.seq.replace(alleles.exp, "G", "R")
	}
	alleles.exp
}

########################################################################################################################

rnb.seq.probe.mismatches <- function(probes.exp, probes.obs) {
	result <- sapply(probes.exp, length) - sapply(probes.obs, length)
	result[sapply(probes.exp, function(x) { length(x) == 1 && is.na(x) })] <- as.integer(NA)
	result[sapply(probes.obs, function(x) { length(x) == 1 && is.na(x) })] <- as.integer(NA)
	i <- which(result < 0L | result > 1L)
	result[i] <- NA
	if (length(i) != 0) {
		warning("Unexpected differences in probe sequence lengths")
	}
	for (i in 1:length(result)) {
		if (is.na(result[i])) {
			next
		}
		x <- probes.exp[[i]]
		if (result[i] == 0) {
			result[i] <- sum(x != probes.obs[[i]])
		} else { # result[i] == 1
			result[i] <- min(sum(x[-length(x)] != probes.obs[[i]]), sum(x[-1] != probes.obs[[i]]))
		}
	}
	result
}

########################################################################################################################

rnb.seq.guess.strand.allele <- function(alleles.exp.pos, alleles.exp.neg, alleles.obs) {
	mismatches.pos <- rnb.seq.probe.mismatches(alleles.exp.pos, alleles.obs)
	mismatches.neg <- rnb.seq.probe.mismatches(alleles.exp.neg, alleles.obs)
	sel.strand <- as.factor(1L + (mismatches.neg < mismatches.pos))
	levels(sel.strand) <- c("+", "-", "*")
	sel.strand[mismatches.pos == mismatches.neg] <- "*"
	mismatches.sel <- ifelse(sel.strand == "+", mismatches.pos, mismatches.neg)
	return(data.frame("Strand" = sel.strand, "Mismatches" = mismatches.sel))
}

########################################################################################################################

rnb.seq.guess.strands <- function(chrom.sequence, loci, pr.design, alleles.A, alleles.B = NULL) {
	alleles.exp.pos <- rnb.seq.get.expected.sequence(chrom.sequence, loci, pr.design, "+")
	alleles.exp.neg <- rnb.seq.get.expected.sequence(chrom.sequence, loci, pr.design, "-")
	alleles.char <- strsplit(alleles.A, NULL)
	if (pr.design == "I") {
		result <- rnb.seq.guess.strand.allele(alleles.exp.pos$A, alleles.exp.neg$A, alleles.char)
		alleles.char <- strsplit(alleles.B, NULL)
		resultB <- rnb.seq.guess.strand.allele(alleles.exp.pos$B, alleles.exp.neg$B, alleles.char)
		i <- which(result$Strand != resultB$Strand)
		result[i, "Strand"] <- "*"
		result <- data.frame(
			"Guessed Strand" = result$Strand,
			"Mismatches A" = result$Mismatches,
			"Mismatches B" = resultB$Mismatches, check.names = FALSE)
	} else { # pr.design == "II"
		result <- rnb.seq.guess.strand.allele(alleles.exp.pos, alleles.exp.neg, alleles.char)
		result <- data.frame(
			"Guessed Strand" = result$Strand,
			"Mismatches A" = result$Mismatches,
			"Mismatches B" = 0L, check.names = FALSE)
	}
	result
}

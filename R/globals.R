########################################################################################################################
## globals.R
## created: 2012-08-18
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Global variables used by functions related to the update of probe annotations.
########################################################################################################################

########################################################################################################################
## Global variables for the targeted genome assembly

## Environment used to store the global variables:
##   assembly, DIR.PACKAGE, GENOME, CHROMOSOMES, snps, regions, sites, mappings, ...
.globals <- new.env()

########################################################################################################################
## Predefined constants used in RnBeads

## Dinucleotide patterns to be annotated
NUCLEOTIDE.PATTERNS <- c("CG")
names(NUCLEOTIDE.PATTERNS) <- sapply(strsplit(NUCLEOTIDE.PATTERNS, ""), paste, collapse = "p")

## Length, in base pairs, of a neighborhood of targeted CpG dinucleotide; used in calculating GC content and CpG density
LENGTH.NEIGHBORHOOD <- 100L

## Length, in base pairs, of tiling regions
LENGTH.TILING <- 5000L

## Length, in base pairs, of a gene promoter region upstream of the TSS
LENGTH.PROMOTER.UPSTREAM <- 1500L

## Length, in base pairs, of a gene promoter region downstream of the TSS
LENGTH.PROMOTER.DOWNSTREAM <- 500L

## Maximum distance, in base pairs, from a CpG island that defines a CGI shore
LENGTH.CGI.SHORE <- 2000L

## Maximum distance, in base pairs, from a CpG island shore that defines a CGI shelf
LENGTH.CGI.SHELF <- 2000L

## Annotation descriptions; included as attributes to the annotation structures
ANNOT.DESCRIPTIONS <- c(
	"CpG" = "CpG dinucleotides",
	"probes450" = "HumanMethylation450 BeadChip probes",
	"tiling" = paste("Genome tiling regions of length", LENGTH.TILING),
	"genes" = "Ensembl genes",
	"promoters" = "Promoter regions of Ensembl genes",
	"cpgislands" = "CpG island track of the UCSC Genome browser")

########################################################################################################################
## dbSNP

## Base FTP location of dbSNP
DBSNP.FTP.BASE <- "ftp://ftp.ncbi.nih.gov/snp/organisms/"

REFERENCE2ASSEMBLY <- c(
	"GRCh38" = "hg38",
	"GRCh37.p13" = "hg19",
	"GCF_000001635.21" = "mm10", # GRCm38.p1
	"GCF_000001635.22" = "mm10", # GRCm38.p2
	"GCF_000001895.4" = "rn5",
	"GCF_000001895.5" = "rn6",
	"GCF_000002035.4" = "zv9")

## Maximum value for a major allele frequency to consider
MAJOR.ALLELE.FREQUENCY <- 0.95

########################################################################################################################
## UCSC Genome Browser

## Base FTP location of files from the UCSC Genome Browser
UCSC.FTP.BASE <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/"

########################################################################################################################
## Gene Expression Omnibus

## Base FPT directory of GEO's records dedicated to HumanMethylation27 and HumanMethylation450
GEO.FTP.BASE <- "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/platforms/"

########################################################################################################################
## Infinium 27k, 450k, and EPIC

## Column names (assigned by this script) of the table on control probes downloaded from GEO and Illumina's site
INF.CONTROL.PROBE.TABLE.COLUMNS <- c("ID", "Target", "Color", "Description", "AVG")

## Control probe colors associated with the evaluation of the Red channel
CONTROL.COLORS.GREEN <- c("Black", "Blue", "Cyan", "Green", "Lime", "LimeGreen", "SkyBlue")

## Control probe colors associated with the evaluation of the Red channel
CONTROL.COLORS.RED <- c("Gold", "Orange", "Purple", "Red", "Tomato", "Yellow")

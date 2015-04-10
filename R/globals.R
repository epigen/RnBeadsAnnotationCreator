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
##   assembly, DIR.PACKAGE, CHROMOSOMES, snps, regions, sites, mappings
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
	"GCF_000001635.21" = "mm10", # Genome Reference Consortium Mouse Build 38 patch release 1 (GRCm38.p1) 
	"GCF_000001895.4" = "rn5")

## Maximum value for a major allele frequency to consider
MAJOR.ALLELE.FREQUENCY <- 0.95
	
########################################################################################################################
## UCSC Genome Browser

## Base FTP location of files from the UCSC Genome Browser
UCSC.FTP.BASE <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/"

## Text file containing the CpG Island track in the UCSC Genome Browser
#UCSC.FTP.CGIS <- paste0(UCSC.FTP.BASE, assembly, "/database/cpgIslandExt.txt.gz")

## Text file containing the definition of chromosome bands
#UCSC.FTP.BANDS <- paste0(UCSC.FTP.BASE, assembly, "/database/cytoBand.txt.gz")
#rm(UCSC.FTP.BASE)

########################################################################################################################
## Gene Expression Omnibus

## Base FPT directory of GEO's records dedicated to HumanMethylation27K and HumanMethylation450K
GEO.FTP.BASE <- "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/platforms/"

########################################################################################################################
## Infinium 27k and 450k

## Column names (assigned by this script) of the table on control probes downloaded from GEO
GEO.CONTROL.PROBE.TABLE.COLUMNS <- c("ID", "Target", "Color", "Description", "AVG")

## Control probe colors associated with the evaluation of the Red channel
CONTROL.COLORS.GREEN <- c("Black", "Blue", "Cyan", "Green", "Lime", "LimeGreen", "SkyBlue")

## Control probe colors associated with the evaluation of the Red channel
CONTROL.COLORS.RED <- c("Gold", "Orange", "Purple", "Red", "Tomato", "Yellow")

########################################################################################################################

## Location of the control probe annotation table compiled by Pavlo Lutsik
COMPILED.TABLE.FILE <- "tables/control.meta.data.csv"

## Column names in the CSV file compiled by Pavlo Lutsik
COMPILED.TABLE.COLUMNS <- c("Address", "Purpose", "Description", "Evaluate.Green", "Evaluate.Red", "Expected.Intensity",
	"Control.Type", "TargetID", "Index")


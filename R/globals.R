########################################################################################################################
## globals.R
## created: 2012-08-18
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Global variables used by functions related to the update of probe annotations.
########################################################################################################################

########################################################################################################################
## Global variables for the targeted genome assembly

## Environment containing the global variables:
##   assembly, DIR.PACKAGE
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

########################################################################################################################
## dbSNP

## Base FTP location of dbSNP
DBSNP.FTP.BASE <- "ftp://ftp.ncbi.nih.gov/snp/organisms/"

## FTP location of BED files for human in dbSNP
#DBSNP.FTP.VCF <- list(
#	"hg19" = paste0(DBSNP.FTP.BASE, "human_9606/VCF/", "00-All.vcf.gz"),
#	"mm10" = paste0(DBSNP.FTP.BASE, "mouse_10090/VCF/", "vcf_chr_", c(1:19, "X", "Y"), ".vcf.gz"),
#	"rn5" = paste0(DBSNP.FTP.BASE, "rat_10116/VCF/", "vcf_chr_", c(1:20, "X"), ".vcf.gz"))
#rm(DBSNP.FTP.BASE)

########################################################################################################################
## UCSC Genome Browser

## Base FTP location of files from the UCSC Genome Browser
UCSC.FTP.BASE <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/"

## Text file containing the CpG Island track in the UCSC Genome Browser
#UCSC.FTP.CGIS <- paste0(UCSC.FTP.BASE, assembly, "/database/cpgIslandExt.txt.gz")

## Text file containing the definition of chromosome bands
#UCSC.FTP.BANDS <- paste0(UCSC.FTP.BASE, assembly, "/database/cytoBand.txt.gz")
#rm(UCSC.FTP.BASE)

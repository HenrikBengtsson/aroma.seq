############################################################################
# DESCRIPTION:
# This script shows how to (1) align single-end DNAseq reads in FASTQ files
# to the human genome (FASTA file), (2) count the aligned reads (BAM files)
# in uniformely distributed 50kb bins (UGP file), (3) normalize the counts
# for amount of GC-content in each bins (UNC file), and (4) finally
# segment the normalized DNAseq total copy-number counts using
# Circular Binary Segmentation (CBS) and (5) generate an interactive
# Chromosome Explorer report viewable in the browser.
#
# REQUIREMENTS:
# fastqData/
#  AlbertsonD_2012-SCC/
#   HomoSapiens/
#    <sample>_<barcode>_L[0-9]{3}_R[12]_[0-9]{3}.fastq [private data]
#
# annotationData/
#  chipTypes/
#   GenericHuman/
#    GenericHuman,50kb,HB20090503.ugp [1]
#    GenericHuman,50kb,HB20121021.unc [1]
#  organisms/
#   HomoSapiens/
#    human_g1k_v37.fasta [2]
#
# REFERENCES:
# [1] http://aroma-project.org/data/annotationData/chipTypes/GenericHuman/
# [2] ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/
#                                                   human_g1k_v37.fasta.gz
############################################################################
library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
fullTest <- fullTest && isPackageInstalled("BatchJobs")
if (fullTest) {
# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "YeastTest"
organism <- "SaccharomycesCerevisiae"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build Bowtie2 index set (iff missing)
is <- buildBowtie2IndexSet(fa, verbose=verbose)
print(is)

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(fqs, SamReadGroup(LB="MPS-034"))

alg <- Bowtie2Alignment(fqs, indexSet=is)
print(alg)

bams <- process(alg, verbose=verbose)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate target UGP on the fly
# NOTE: This is a rather ad hoc solution in order to be able
#       to tie into the aroma.cn framework. /HB 2013-11-12
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chrLengths <- getSeqLengths(fa);
by <- 25e3
chrBins <- lapply(chrLengths, FUN=function(n) seq(from=1L, to=n, by=by))
nbrOfLoci <- sum(sapply(chrBins, FUN=length))

byTag <- sprintf("%gkb", by/1e3)
filename <- sprintf("%s.ugp", paste(c(organism, byTag), collapse=","))
ugp <- AromaUgpFile$allocate(filename=filename, path=tempfile(), nbrOfRows=nbrOfLoci, platform="HT-Seq", chipType=organism, overwrite=TRUE);
offset <- 0L
for (chr in seq_along(chrBins)) {
  x <- chrBins[[chr]]
  idxs <- offset + seq_along(x)
  ugp[idxs,1] <- chr
  ugp[idxs,2] <- x
  offset <- offset + length(x)
}
print(ugp)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads per bin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- TotalCnBinnedCounting(bams, targetUgp=ugp)
print(bc)

counts <- process(bc, verbose=verbose)
verbose && print(verbose, counts)


} # if (fullTest)

############################################################################
# HISTORY:
# 2013-11-11
# o Created from BinnedGcNormalization test script.
############################################################################

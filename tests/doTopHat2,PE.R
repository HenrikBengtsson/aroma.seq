library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "samtools")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "YeastTest"
organism <- "SaccharomycesCerevisiae"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# FASTQ data
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doTopHat2(fqs, reference=fa, transcripts=NULL, verbose=TRUE)
print(bams)

# Specifying additional parameters (identical to the defaults)
bams2 <- doTopHat2(fqs, reference=fa, transcripts=NULL, mateInnerDist=50, mateStdDev=20, tags=c("*", "50_20"), verbose=TRUE)
print(bams2)

# Assert that we get the same as the default settings
# Comparison of BAM files must be done at the SAM level here,
# because there are additional things encoded in the BAM files.
sams <- convertToSam(bams)
sams2 <- convertToSam(bams2)
samsC <- getChecksumFileSet(sams)
sams2C <- getChecksumFileSet(sams2)
print(samsC)
print(sams2C)
#  stopifnot(equals(samsC, sams2C))

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-11-02
# o Created.
############################################################################

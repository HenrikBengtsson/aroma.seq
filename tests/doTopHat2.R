library("aroma.seq")
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
bams <- doTopHat2(fqs, reference=fa, verbose=TRUE)
print(bams)

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-11-02
# o Created.
############################################################################

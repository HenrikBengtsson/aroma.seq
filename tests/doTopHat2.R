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

# Setup (writable) local data directory structure
setupExampleData()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
path <- file.path("annotationData", "organisms", organism)
fas <- FastaReferenceSet$byPath(path)
fa <- getFile(fas, 1)  # Presuming there is only one reference fasta file
print(fa)

# FASTQ data
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path, paired=TRUE)
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

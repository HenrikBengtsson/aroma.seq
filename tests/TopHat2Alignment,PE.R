library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

dataSet <- "YeastTest"
organism <- "Saccharomyces_cerevisiae"

# Setup (writable) local data directory structure
setupExampleData()


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
is <- buildBowtie2IndexSet(fa, verbose=TRUE)  # is = 'index set'
print(is)

# Align input reads using TopHat
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
print(ta)

fullTest <- fullTest && isCapableOf(aroma.seq, "samtools")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {
process(ta, verbose=-100)

bams <- getOutputDataSet(ta)
print(bams)

# Sanity checks
stopifnot(length(bams) == length(fqs))
} # if (fullTest)

} # if (fullTest)

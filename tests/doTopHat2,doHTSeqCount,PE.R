library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
fullTest <- fullTest && isCapableOf(aroma.seq, "htseq")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "YeastTest"
organism <- "Saccharomyces_cerevisiae"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

gtf <- GtfDataFile$byOrganism(organism)
print(gtf)

# FASTQ data
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)

# Make sure to group FASTQ files per sample
fqs <- setFullNamesTranslator(fqs, function(names, ...) {
  names <- gsub("_(1|R1)$", "", names)
  gsub("_", ",", names, fixed=TRUE)
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doTopHat2(fqs, reference=fa, transcripts=gtf, groupBy="name", verbose=TRUE)
print(bams)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HTSeqCount
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
counts <- doHTSeqCount(bams, transcripts=gtf, verbose=-100)
print(counts)

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-01-24
# o Created.
############################################################################

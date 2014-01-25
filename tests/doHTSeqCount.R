library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
fullTest <- fullTest && isCapableOf(aroma.seq, "htseq")
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

gtf <- GtfDataFile$byOrganism(organism)
print(gtf)

# FASTQ data
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)

# Set fullnames translator, making SRR + 5 digits the name
# and the rest tags, just as an example
fqs <- setFullNamesTranslator(fqs, function(names, ...) {
  # Drop any stray "R1" suffix
  names <- gsub("_(1|R1)$", "", names)
  # Tagify
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
counts <- doHTSeqCount(bams, transcripts=gtf, verbose=TRUE)
print(counts)

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-01-24
# o Created.
############################################################################

library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

dataSet <- "YeastTest"
organism <- "SaccharomycesCerevisiae"

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


# Set fullnames translator, making SRR + 5 digits the name
# and the rest tags, just as an example
fqs <- setFullNamesTranslator(fqs, function(names, ...) {
  # Drop any stray "R1" suffix
  names <- gsub("_(1|R1)$", "", names)
  # Tagify
  gsub("_", ",", names, fixed=TRUE)
})
print(getFullNames(fqs))

ta <- TopHat2Alignment(dataSet=fqs, indexSet=is, groupBy="name")
print(ta)

# Assert that for each group a unique set of files are identified
groups <- getGroups(ta)
lapply(groups, FUN=function(idxs) stopifnot(!anyDuplicated(idxs)))


fullTest <- fullTest && isCapableOf(aroma.seq, "samtools")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {
process(ta, verbose=-100)

bams <- getOutputDataSet(ta)
print(bams)

} # if (fullTest)

} # if (fullTest)

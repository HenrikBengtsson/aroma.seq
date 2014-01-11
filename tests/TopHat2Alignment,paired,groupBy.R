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
  # Also, drop the "R1" suffix
  pattern <- "^(SRR[0-9]{5})([0-9])_(chr[0-9]+)_(1|R1)$";
  gsub(pattern, "\\1,sub\\2", names)
})
print(getFullNames(fqs))

ta <- TopHat2Alignment(dataSet=fqs, indexSet=is, groupBy="name")
print(ta)

fullTest <- fullTest && isCapableOf(aroma.seq, "samtools")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {
process(ta, verbose=TRUE)

bams <- getOutputDataSet(ta)
print(bams)

} # if (fullTest)

} # if (fullTest)

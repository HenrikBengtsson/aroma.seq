library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
## Disable test because the GTF and FASTA files are incompatible.
if (fullTest) {

dataset <- "YeastTest"
organism <- "Saccharomyces_cerevisiae"

# Setup (writable) local data directory structure
setupExampleData()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# Bowtie2 index set
is <- buildBowtie2IndexSet(fa, verbose=TRUE)  # is = 'index set'
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FASTQ data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataset, organism=organism, paired=TRUE)
print(fqs)

# Set fullnames translator, making SRR + 5 digits the name
# and the rest tags, just as an example
fqs <- setFullNamesTranslator(fqs, function(names, ...) {
  # Drop any stray "R1" suffix
  names <- gsub("_(1|R1)$", "", names)
  # Tagify
  gsub("_", ",", names, fixed=TRUE)
})
print(getFullNames(fqs))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ta <- TopHat2Alignment(dataSet=fqs, groupBy="name", indexSet=is)
print(ta)

# Assert that for each group a unique set of files are identified
groups <- getGroups(ta)
lapply(groups, FUN=function(idxs) stopifnot(!anyDuplicated(idxs)))

fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {
  bams <- process(ta, verbose=-100)
  print(bams)
} # if (fullTest)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat alignment with transcriptome model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gtf <- GtfDataFile$byOrganism(organism)
print(gtf)
ta <- TopHat2Alignment(dataSet=fqs, groupBy="name", indexSet=is, transcripts=gtf)
print(ta)

if (fullTest) {
  bams <- process(ta, verbose=-100)
  print(bams)
  # Assert proper name and organism inference
  stopifnot(getName(bams) == dataset)
  stopifnot(getOrganism(bams) == organism)
  stopifnot(all(getNames(bams) == names(groups)))
} # if (fullTest)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HT-Seq counting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fullTest <- fullTest && isCapableOf(aroma.seq, "htseq")
if (fullTest) {
  htc <- HTSeqCounting(bams, transcripts=gtf)
  print(htc)

  counts <- process(htc, verbose=-100)
  print(counts)

  # Assert proper name and organism inference
  stopifnot(getName(counts) == dataset)
  stopifnot(getOrganism(counts) == organism)
  stopifnot(all(getNames(counts) == getNames(bams)))

  if (isPackageInstalled("edgeR")) {
    dge <- extractDGEList(counts)
    print(dge)
  }
} # if (fullTest)

} # if (fullTest)

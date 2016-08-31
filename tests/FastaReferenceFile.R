library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "Lambda_phage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build FASTA index file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fai <- buildIndex(fa, verbose=-10)
print(fai)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (fullTest && isCapableOf(aroma.seq, "bwa")) {
  is <- buildBwaIndexSet(fa, verbose=-10)
  print(is)

  # Assert that both ways to build the BWA index generates
  # the exact same set of indices.
  isA <- buildBwaIndexSet(fa, method="is", verbose=-10)
  print(isA)
  isB <- buildBwaIndexSet(fa, method="bwtsw", verbose=-10)
  print(isB)
  checksumA <- sapply(isA, FUN=getChecksum)
  checksumB <- sapply(isB, FUN=getChecksum)
  stopifnot(identical(checksumA, checksumB))
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build Bowtie2 index file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (fullTest && isCapableOf(aroma.seq, "bowtie2")) {
  is <- buildBowtie2IndexSet(fa, verbose=-10)
  print(is)
}



organism <- "Saccharomyces_cerevisiae"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Checksums per sequence
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- getSeqChecksums(fa, verbose=-10)
print(cs)




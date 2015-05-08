library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {


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
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, verbose=-10)
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(fqs, indexSet=is, n=2, q=40)
print(alg)

bams <- process(alg, verbose=-20)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (isCapableOf(aroma.seq, "picard")) {
  # Without IGNORE="MISSING_READ_GROUP" below,
  # we get error 'Read groups is empty'
  bam <- bams[[1]]
  validate(bam, IGNORE="MISSING_READ_GROUP")
  validate(bam, onError="warning")
}


} # if (fullTest)


############################################################################
# HISTORY:
# 2014-04-13
# o Created from corresponding SE test.
############################################################################

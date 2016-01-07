library("aroma.seq")
setOption("R.filesets/parallel", "BiocParallel")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isPackageInstalled("BatchJobs")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

library("future")
strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (require(pkg <- "async", character.only=TRUE)) {
  backend("local")
  strategies <- c(strategies, "batchjobs")
}
setOption("R.filesets/parallel", "future")

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
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (strategy in strategies) {
  plan(strategy)
  print(plan())
  bams <- doBowtie2(fqs, reference=fa, tags=c("*", strategy), verbose=-20)
  print(bams)
}

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-26
# o Created from doBowtie2.R.
############################################################################
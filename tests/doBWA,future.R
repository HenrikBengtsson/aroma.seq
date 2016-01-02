library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {

library("future")
strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (require("async")) {
  strategies <- c(strategies, "batchjobs")
  async::backend("local")
}
setOption("R.filesets/parallel", "future")


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
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (strategy in strategies) {
  plan(strategy)
  print(plan())
  bams <- doBWA(fqs, reference=fa, tags=c("*", strategy), verbose=-20)
  print(bams)
}

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-31
# o Created from doBowtie2,parallel.R.
############################################################################

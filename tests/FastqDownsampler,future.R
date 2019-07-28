library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && (Sys.getenv("_R_CHECK_BUGGY_") != "")
fullTest <- fullTest && isPackageInstalled("ShortRead")
if (fullTest) {

library("future")
strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (require(pkg <- "future.BatchJobs", character.only=TRUE)) {
  strategies <- c(strategies, "batchjobs_local")
}

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "Lambda_phage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (strategy in strategies) {
  plan(strategy)
  print(plan())

  ds <- FastqDownsampler(fqs, subset=25, tags=c("*", strategy))
  print(ds)
  
  fqsS <- process(ds, verbose=-10)
  print(fqsS)

  # Sanity checks
  stopifnot(identical(getFullNames(fqsS), getFullNames(fqs)))
  fqS <- fqsS[[1]]
  stopifnot(nbrOfSeqs(fqS) == getSampleSize(ds))
}

} # if (fullTest)

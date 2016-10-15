library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

library("future")
strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (require(pkg <- "future.BatchJobs", character.only=TRUE)) {
  strategies <- c(strategies, "batchjobs_local")
}

setupExampleData()
dataSet <- "TopHat-example"
organism <- "Lambda_phage"
fa <- FastaReferenceFile$byOrganism(organism)
fqs <- FastqDataSet$byName(dataSet, organism=organism)
bams <- doBowtie2(fqs, reference=fa, verbose=-20)
print(bams)

bams <- setFullNamesTranslator(bams, function(names, ...) {
  sprintf("SampleA,%s", names)
})
print(getFullNames(bams))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up BamMerger:s
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (strategy in strategies) {
  plan(strategy)
  print(plan())

  bm <- BamMerger(bams, groupBy="name", tags=c("*", strategy))
  print(bm)
  groups <- getGroups(bm)
  print(groups)

  # Merging
  bamsM <- process(bm, verbose=-20)
  print(bamsM)
}


} # if (fullTest)

library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism)
print(fqs)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FastQC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rep <- FastQCReporter(fqs)
print(rep)

if (fullTest && isCapableOf(aroma.seq, "fastqc")) {
  res <- process(rep, verbose=-10)
  print(res)
}

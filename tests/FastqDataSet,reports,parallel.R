library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && require("qrqc")

# Setup (writable) local data directory structure
setOption(aromaSettings, "devel/parallel", "BiocParallel::BatchJobs")

setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism)
print(fqs)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate QC report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (fullTest) {
  reports <- report(fqs, verbose=-10)
  print(reports)
}

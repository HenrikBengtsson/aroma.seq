library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "BiocParallel")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && require("qrqc")

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate QC report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (fullTest) {
  html <- report(fqs, verbose=-10)
  print(html)
}

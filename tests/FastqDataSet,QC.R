library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && require("qrqc")

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path)
print(fqs)
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  print(fq)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate QC report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (fullTest) {
  pdfs <- report(fqs, verbose=-10)
  print(pdfs)
}

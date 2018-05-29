library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

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
# FastQC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rep <- FastQCReporter(fqs, groupBy="name")
print(rep)

if (fullTest && isCapableOf(aroma.seq, "fastqc")) {
  res <- process(rep, verbose=-10)
  print(res)
}

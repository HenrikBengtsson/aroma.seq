library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

# Setup (writable) local data directory structure
setupExampleData()


dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism)
print(fqs)
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  print(fq)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gzip data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataSetZ <- sprintf("%s,gz", dataSet);
pathZ <- file.path("fastqData", dataSetZ, organism)
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  pathnameZ <- file.path(pathZ, sprintf("%s.gz", getFilename(fq)))
  gzip(getPathname(fq), pathnameZ, skip=TRUE, remove=FALSE)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup gzip'ed FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqsZ <- FastqDataSet$byName(dataSet, tags="gz", organism=organism)
print(fqsZ)
for (ii in seq_along(fqsZ)) {
  fqZ <- getFile(fqsZ, ii)
  print(fqZ)
}

# Assert that file extensions are properly dropped
stopifnot(all(getFullNames(fqsZ) == getFullNames(fqs)))

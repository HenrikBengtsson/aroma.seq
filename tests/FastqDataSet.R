library("aroma.seq")
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
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path)
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
  if (!isFile(pathnameZ)) gzip(getPathname(fq), pathnameZ, remove=FALSE)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup gzip'ed FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqsZ <- FastqDataSet$byPath(pathZ)
print(fqsZ)
for (ii in seq_along(fqsZ)) {
  fqZ <- getFile(fqsZ, ii)
  print(fqZ)
}

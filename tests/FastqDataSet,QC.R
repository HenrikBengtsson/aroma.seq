library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && require("qrqc")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup (writable) local data directory structure
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathD <- system.file("exData", package="aroma.seq")
for (dir in c("fastqData")) {
  copyDirectory(file.path(pathD, dir), to=dir, overwrite=FALSE)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData", "TopHat-example", "LambdaPhage")
fqs <- FastqDataSet$byPath(path, pattern="[.](fq|fastq)$")
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

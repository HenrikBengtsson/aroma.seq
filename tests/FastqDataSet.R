library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

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
# Gzip data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathGZ <- file.path("fastqData", "TopHat-example,gz", "LambdaPhage")
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  pathnameGZ <- file.path(pathGZ, sprintf("%s.gz", getFilename(fq)))
  if (!isFile(pathnameGZ)) {
    fqGZ <- gzip(getPathname(fq), pathnameGZ)
  }
  print(fqGZ)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup gzip'ed FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- pathGZ
fqs <- FastqDataSet$byPath(path, pattern="[.]gz$")
print(fqs)
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  print(fq)
}

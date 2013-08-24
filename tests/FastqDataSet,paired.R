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
# Setup paired-end FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData", "TopHat-example", "LambdaPhage")
fqs <- FastqDataSet$byPath(path, pattern="_1[.](fq|fastq)$", paired=TRUE)
print(fqs)
pairs <- getFilePairs(fqs)
print(pairs)
print(pairs[,1])
print(pairs[,2])


# Locate mate pair
r1 <- getFile(fqs, 1)
print(r1)
r2 <- getMateFile(r1)
print(r2)
r1b <- getMateFile(r2)
stopifnot(equals(r1b, r1))


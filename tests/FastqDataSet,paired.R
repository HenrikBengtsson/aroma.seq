library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

# Setup (writable) local data directory structure
setupExampleData()


dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup paired-end FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
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
print(r1b)
stopifnot(identical(getPathname(r1b), getPathname(r1)))

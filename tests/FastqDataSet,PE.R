library("aroma.seq")

# Setup (writable) local data directory structure
setupExampleData()


dataSet <- "TopHat-example"
organism <- "Lambda_phage"

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
r1 <- fqs[[1]]
print(r1)
r2 <- getMateFile(r1)
print(r2)
r1b <- getMateFile(r2)
print(r1b)
stopifnot(identical(getPathname(r1b), getPathname(r1)))

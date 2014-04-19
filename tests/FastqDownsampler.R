library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && (Sys.getenv("_R_CHECK_BUGGY_") != "")
fullTest <- fullTest && isPackageInstalled("ShortRead")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)

n <- sapply(fqs, FUN=nbrOfSeqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample to fixed count
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- FastqDownsampler(fqs, subset=25)
print(ds)

fqsS <- process(ds, verbose=-10)
print(fqsS)

# Sanity checks
stopifnot(identical(getFullNames(fqsS), getFullNames(fqs)))
nS <- sapply(fqsS, FUN=nbrOfSeqs)
stopifnot(all(nS == 25L))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample to fraction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- FastqDownsampler(fqs, subset=0.10)
print(ds)

fqsS <- process(ds, verbose=-20)
print(fqsS)

# Sanity checks
stopifnot(identical(getFullNames(fqsS), getFullNames(fqs)))
nS <- sapply(fqsS, FUN=nbrOfSeqs)
stopifnot(all(nS == 0.10*n))

} # if (fullTest)

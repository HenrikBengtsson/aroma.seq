library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
## fullTest <- fullTest && (Sys.getenv("_R_CHECK_BUGGY_") != "")
fullTest <- fullTest && isPackageInstalled("ShortRead")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "Lambda_phage"

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

md5S0 <- sapply(fqsS, FUN=getChecksum)
print(md5S0)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reproducibility
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqsS <- doDownsample(fqs, subset=0.10, seed=0xBEEF, tags=c("*", "seed=0xBEEF"))
print(fqsS)
md5S1 <- sapply(fqsS, FUN=getChecksum)
print(md5S1)

fqsS <- doDownsample(fqs, subset=0.10, seed=0xBEEF, tags=c("*", "seed=0xBEEF", "r2"))
print(fqsS)
md5S2 <- sapply(fqsS, FUN=getChecksum)
print(md5S2)

# Sanity check
stopifnot(all(md5S2 == md5S1))

} # if (fullTest)

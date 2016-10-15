library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && (Sys.getenv("_R_CHECK_BUGGY_") != "")
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- FastqDownsampler(fqs, subset=25, tags=c("*", "parallel"))
print(ds)

fqsS <- process(ds, verbose=-10)
print(fqsS)

# Sanity checks
stopifnot(identical(getFullNames(fqsS), getFullNames(fqs)))
fqS <- fqsS[[1]]
stopifnot(nbrOfSeqs(fqS) == getSampleSize(ds))

} # if (fullTest)

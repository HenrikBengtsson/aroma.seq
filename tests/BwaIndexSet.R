library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {

## Setup (writable) local data directory structure
setupExampleData()

dataset <- "YeastTest"
organism <- "Saccharomyces_cerevisiae"

## Setup FASTA reference file
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

## Build index set
is <- buildBwaIndexSet(fa, verbose=-10)
print(is)

stopifnot(isCompatibleWith(is, fa))
stopifnot(isCompatibleWith(fa, is))

} # if (fullTest)

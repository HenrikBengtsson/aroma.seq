library("aroma.seq")
equals <- R.oo::equals # In case testthat::equals() is attached

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

setupExampleData()

organism <- "Saccharomyces_cerevisiae"
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

fa2 <- FastaReferenceFile$byOrganism(organism, prefix="SC_chr1-2")
print(fa2)

# Sanity check
stopifnot(equals(fa2, fa))

is <- buildBowtie2IndexSet(fa, verbose=TRUE)
print(is)

fa3 <- getFastaReferenceFile(is)
print(fa3)

# Sanity check
stopifnot(equals(fa3, fa))


} # if (fullTest)

library("aroma.seq")

if (isCapableOf(aroma.seq, "gatk")) {
  bin <- findGATK()
  print(bin)


}

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "gatk")
fullTest <- fullTest && isCapableOf(aroma.seq, "picard")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {


# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

fai <- buildIndex(fa)
print(fai)

dict <- buildDictionary(fa)
print(dict)

fqs <- FastqDataSet$byName(dataSet, organism=organism)
print(fqs)

bams <- doBowtie2(fqs, reference=fa, verbose=TRUE)
print(bams)

pathnameFA <- getPathname(fa)
pathnameBAMs <- getPathnames(bams)

tryCatch({
  res <- gatk(analysisType="CountReads", pathnameR=pathnameFA,
              pathnameI=pathnameBAMs, verbose=TRUE)
  print(res)
}, error = function(ex) {
  print(ex)
})

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-04-13
# o Created.
############################################################################

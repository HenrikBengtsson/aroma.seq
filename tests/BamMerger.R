library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

setupExampleData()
dataSet <- "TopHat-example"
organism <- "LambdaPhage"
fa <- FastaReferenceFile$byOrganism(organism)
fqs <- FastqDataSet$byName(dataSet, organism=organism)
bams <- doBowtie2(fqs, reference=fa, verbose=-20)
print(bams)

bams <- setFullNamesTranslator(bams, function(names, ...) {
  sprintf("SampleA,%s", names)
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up BamMerger:s
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bm <- BamMerger(bams, groupBy="name")
print(bm)
groups <- getGroups(bm)
print(groups)

bm2 <- BamMerger(bams, groupBy=groups)
print(bm2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merging
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bamsM <- process(bm, verbose=-20)
print(bamsM)


} # if (fullTest)

library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

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
print(getFullNames(bams))

n <- sapply(bams, FUN=nbrOfReads)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample to fixed counts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- BamDownsampler(bams, subset=25)
print(ds)

bamsS <- process(ds, verbose=-20)
print(bamsS)

# Sanity checks
stopifnot(identical(getFullNames(bamsS), getFullNames(bams)))
nS <- sapply(bamsS, FUN=nbrOfReads)
stopifnot(all(nS == 25L))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Downsample to fraction
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- BamDownsampler(bams, subset=0.10)
print(ds)

bamsS <- process(ds, verbose=-20)
print(bamsS)

# Sanity checks
stopifnot(identical(getFullNames(bamsS), getFullNames(bams)))
nS <- sapply(bamsS, FUN=nbrOfReads)
stopifnot(all(nS == 0.10*n))

} # if (fullTest)

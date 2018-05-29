library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

setupExampleData()
dataSet <- "TopHat-example"
organism <- "Lambda_phage"
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

md5S0 <- sapply(bamsS, FUN=getChecksum)
print(md5S0)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reproducibility
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bamsS <- doDownsample(bams, subset=0.10, seed=0xBEEF, tags=c("*", "seed=0xBEEF"))
print(bamsS)
md5S1 <- sapply(bamsS, FUN=getChecksum)
print(md5S1)

bamsS <- doDownsample(bams, subset=0.10, seed=0xBEEF, tags=c("*", "seed=0xBEEF", "r2"))
print(bamsS)
md5S2 <- sapply(bamsS, FUN=getChecksum)
print(md5S2)

# Sanity check
stopifnot(all(md5S2 == md5S1))

} # if (fullTest)

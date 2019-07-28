library("aroma.seq")

fullTest <- isCapableOf(aroma.seq, "bwa")
if (fullTest) {


# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "Lambda_phage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doBWA(fqs, reference=fa, verbose=-200)
print(bams)

bams2 <- doBWA(bams, reference=fa, verbose=-200)
print(bams2)


bam <- bams[[1]]
print(bam)

bam2 <- bams2[[1]]
print(bam2)


## Assert that the realigned files contains identical alignment
sam <- convertToSam(bam, overwrite = TRUE)
print(sam)
sam2 <- convertToSam(bam2, overwrite = TRUE)
print(sam2)
bfr <- readLines(getPathname(sam))
bfr2 <- readLines(getPathname(sam2))
bfr <- grep("@PG", bfr, value = TRUE, invert = TRUE)
bfr2 <- grep("@PG", bfr2, value = TRUE, invert = TRUE)
stopifnot(all.equal(bfr2, bfr))

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-22
# o Created.
############################################################################

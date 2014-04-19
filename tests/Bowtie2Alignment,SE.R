library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, verbose=-10)
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
alg <- Bowtie2Alignment(fqs, indexSet=is)
print(alg)

bams <- process(alg, verbose=-20)
print(bams)

# Display an example BAM file
for (ii in seq_along(bams)) print(bams[[ii]])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment on gzip'ed FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gzip data set
dataSetZ <- sprintf("%s,gz", dataSet);
pathZ <- file.path("fastqData", dataSetZ, organism)
for (ii in seq_along(fqs)) {
  fq <- fqs[[ii]]
  pathnameZ <- file.path(pathZ, sprintf("%s.gz", getFilename(fq)))
  if (!isFile(pathnameZ)) gzip(getPathname(fq), pathnameZ, remove=FALSE)
}
fqsZ <- FastqDataSet$byName(dataSet, tags="gz", organism=organism)

# Alignment
algZ <- Bowtie2Alignment(fqsZ, indexSet=is)
print(algZ)

bamsZ <- process(algZ, verbose=-20)
print(bamsZ)

# Results should be identical with and without gzip'ed FASTQ files
stopifnot(length(bamsZ) == length(bams))
stopifnot(identical(getFullNames(bamsZ), getFullNames(bams)))
for (ii in seq_along(bams)) {
  bam <- bams[[ii]]
  bamZ <- bamsZ[[ii]]
  stopifnot(getChecksum(bamZ) == getChecksum(bam))
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (isCapableOf(aroma.seq, "picard")) {
  bam <- bams[[1]]
  # Without IGNORE="MISSING_READ_GROUP" we get an
  # error on 'Read groups is empty'
  validate(bam, IGNORE="MISSING_READ_GROUP")
  validate(bam, onError="warning")
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove duplicated reads using Picard
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (isCapableOf(aroma.seq, "picard")) {
  dr <- PicardDuplicateRemoval(bams)
  print(dr)

  bamsU <- process(dr, verbose=-20)
  print(bamsU)
}


} # if (fullTest)


############################################################################
# HISTORY:
# 2012-09-27
# o Created.
############################################################################

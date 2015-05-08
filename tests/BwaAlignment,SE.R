library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
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
# Build index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, verbose=-10)
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(fqs, indexSet=is, n=2, q=40)
print(alg)

bams <- process(alg, verbose=-20)
print(bams)

# Display an example BAM file
for (ii in seq_along(bams)) print(bams[[ii]])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A side-effect of BWA is for now that SAM data files are also created
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sams <- SamDataSet$byPath(getPath(bams))
print(sams)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (isCapableOf(aroma.seq, "picard")) {
  # Without IGNORE="MISSING_READ_GROUP" below,
  # we get error 'Read groups is empty'
  sam <- sams[[1]]
  bam <- bams[[1]]
  validate(sam, IGNORE="MISSING_READ_GROUP")
  validate(sam, onError="warning")
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment on gzip'ed FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gzip data set
dataSetZ <- sprintf("%s,gz", dataSet);
pathZ <- file.path("fastqData", dataSetZ, organism);
for (ii in seq_along(fqs)) {
  fq <- fqs[[ii]]
  pathnameZ <- file.path(pathZ, sprintf("%s.gz", getFilename(fq)))
  if (!isFile(pathnameZ)) gzip(getPathname(fq), pathnameZ, remove=FALSE)
}
fqsZ <- FastqDataSet$byName(dataSet, tags="gz", organism=organism, paired=FALSE)

# BWA with BWA 'aln' options '-n 2' and '-q 40'.
algZ <- BwaAlignment(fqsZ, indexSet=is, n=2, q=40)
print(algZ)

bamsZ <- process(algZ, verbose=-20)
print(bamsZ)

# Results should be identical with and without gzip'ed FASTQ files
stopifnot(length(bamsZ) == length(bams))
stopifnot(identical(getFullNames(bamsZ), getFullNames(bams)))
for (ii in seq_along(bams)) {
  bam <- bams[[ii]]
  bamZ <- bamsZ[[ii]]
##  stopifnot(getChecksum(bamZ) == getChecksum(bam))
}

} # if (fullTest)


############################################################################
# HISTORY:
# 2012-09-27
# o Created.
############################################################################

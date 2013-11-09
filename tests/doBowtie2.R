library("aroma.seq")
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
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doBowtie2(fqs, reference=fa, verbose=-20)
print(bams)

# Display individual BAM files
for (ii in seq_along(bams)) {
  bam <- getFile(bams, ii)
  print(bam)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment on gzip'ed FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gzip data set
dataSetZ <- sprintf("%s,gz", dataSet);
pathZ <- file.path("fastqData", dataSetZ, organism);
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  pathnameZ <- file.path(pathZ, sprintf("%s.gz", getFilename(fq)))
  if (!isFile(pathnameZ)) gzip(getPathname(fq), pathnameZ, remove=FALSE)
}
fqsZ <- FastqDataSet$byPath(pathZ)

bamsZ <- doBowtie2(fqsZ, reference=fa, verbose=-20)
print(bamsZ)


# Results should be identical with and without gzip'ed FASTQ files
stopifnot(length(bamsZ) == length(bams))
stopifnot(identical(getFullNames(bamsZ), getFullNames(bams)))
for (ii in seq_along(bams)) {
  bam <- getFile(bams, ii)
  bamZ <- getFile(bamsZ, ii)
  stopifnot(getChecksum(bamZ) == getChecksum(bam))
}


} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-22
# o Created.
############################################################################

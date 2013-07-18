library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup (writable) local data directory structure
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathD <- system.file("exData", package="aroma.seq")
for (dir in c("annotationData", "fastqData")) {
  copyDirectory(file.path(pathD, dir), to=dir, overwrite=FALSE)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "organisms", "LambdaPhage")
fa <- FastaReferenceFile("lambda_virus.fa", path=path)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData", "TopHat-example", "LambdaPhage")
fqs <- FastqDataSet$byPath(path, pattern="[.](fq|fastq)$")
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build Bowtie2 index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, verbose=-10)
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end Bowtie2 alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
alg <- Bowtie2Alignment(fqs, indexSet=is)
print(alg)

bams <- process(alg, verbose=-20)
print(bams)

# Display an example BAM file
for (ii in seq_along(bams)) {
  bam <- getFile(bams, ii)
  print(bam)
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
# Single-end Bowtie2 alignment on gzip'ed FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gzip data set
pathZ <- file.path("fastqData", "TopHat-example,gz", "LambdaPhage")
for (ii in seq_along(fqs)) {
  fq <- getFile(fqs, ii)
  pathnameZ <- file.path(pathZ, sprintf("%s.gz", getFilename(fq)))
  if (!isFile(pathnameZ)) gzip(getPathname(fq), pathnameZ, remove=FALSE)
}
fqsZ <- FastqDataSet$byPath(pathZ, pattern="[.]gz$")

# Bowtie2 alignment
algZ <- Bowtie2Alignment(fqsZ, indexSet=is)
print(algZ)

bamsZ <- process(algZ, verbose=-20)
print(bamsZ)

# Bowtie2 results should be identical with and without gzip'ed FASTQ files
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
# 2012-09-27
# o Created.
############################################################################

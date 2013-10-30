library("aroma.seq")

test <- isCapableOf(aroma.seq, "bowtie2")
test <- test && isCapableOf(aroma.seq, "tophat2")

if (test) {
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Test for non-compatible bowtie2 and tophat2 versions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verT <- attr(findTopHat2(), "version")
verB <- attr(findBowtie2(), "version")
bad <- (verT == "2.0.3" && verB == "2.1.0")
if (bad) {
  throw(sprintf("TopHat2 v%s is known to not work with Bowtie2 v%s.", verT, verB))
}


dataSet <- "YeastTest"
organism <- "SC"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create local unzipped copies of data directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The directory where all example data files are
path0 <- system.file(file.path("exData"), package="aroma.seq", mustWork=TRUE)

# Annotation data
path <- file.path("annotationData", "organisms", organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)
sapply(GenericDataFileSet$byPath(path, pattern="[.]gz$"), FUN=gunzip)


# FASTQ data
path <- file.path("fastqData", dataSet, organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)
sapply(GenericDataFileSet$byPath(path, pattern="[.]gz$"), FUN=gunzip)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
path <- file.path("annotationData", "organisms", organism)
fas <- FastaReferenceSet$byPath(path)
fa <- getFile(fas, 1)  # Presuming there is only one reference fasta file
print(fa)

# FASTQ data
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, verbose=TRUE)  # is = 'index set'
print(is)

# Align input reads using TopHat
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
process(ta, verbose=TRUE)

bams <- getOutputDataSet(ta)
print(bams)

# Sanity checks
stopifnot(length(bams) == length(fqs))

} # if (test)

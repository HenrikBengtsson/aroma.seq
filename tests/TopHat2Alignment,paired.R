library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
if (fullTest) {

dataSet <- "YeastTest"
organism <- "SC"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create local unzipped copies of data directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The directory where all example data files are
path0 <- system.file(file.path("exData"), package="aroma.seq", mustWork=TRUE)

# FASTA data
path <- file.path("annotationData", "organisms", organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)

# FASTQ data
path <- file.path("fastqData", dataSet, organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)


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

} # if (fullTest)

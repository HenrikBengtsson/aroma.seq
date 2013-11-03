library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "samtools")
fullTest <- fullTest && isCapableOf(aroma.seq, "tophat2")
fullTest <- fullTest && isPackageInstalled(aroma.seq, "BatchJobs")
if (fullTest) {

setOption(aromaSettings, "devel/parallel", "BiocParallel::BatchJobs")

dataSet <- "YeastTest"
organism <- "SC"

# Setup (writable) local data directory structure
setupExampleData()


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
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is, tags=c("*", "parallel"))
process(ta, verbose=TRUE)

bams <- getOutputDataSet(ta)
print(bams)

# Sanity checks
stopifnot(length(bams) == length(fqs))

} # if (fullTest)

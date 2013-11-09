library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
fullTest <- fullTest && isPackageInstalled("BatchJobs")
if (fullTest) {

setOption(aromaSettings, "devel/parallel", "BiocParallel::BatchJobs")

# Setup (writable) local data directory structure
setupExampleData()


dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "organisms", organism)
fa <- FastaReferenceFile("lambda_virus.fa", path=path)
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
bams <- doBWA(fqs, reference=fa, tags=c("*", "parallel"), verbose=-20)
print(bams)

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-31
# o Created from doBowtie2,parallel.R.
############################################################################

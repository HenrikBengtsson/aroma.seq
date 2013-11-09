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
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism)
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

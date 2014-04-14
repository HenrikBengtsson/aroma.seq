library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
fullTest <- fullTest && isCapableOf(aroma.seq, "gatk")
fullTest <- fullTest && isCapableOf(aroma.seq, "picard")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataset <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

fai <- buildIndex(fa)
dict <- buildDictionary(fa)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataset, organism=organism)
print(fqs)

# Set read groups ('SM' is required by GATK)
lapply(fqs, FUN=function(fq) {setSamReadGroup(fq, SamReadGroup(
  LB=getFullName(fqs),  # Library
  SM=getFullName(fq)    # Sample
))});


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doBowtie2(fqs, reference=fa, verbose=TRUE)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call GATK CountReads, cf. '(howto) Run the GATK for the first time'
# http://gatkforums.broadinstitute.org/discussion/1209/howto-run-the-gatk-for-the-first-time
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathnameFA <- getPathname(fa)
pathnameBAMs <- getPathnames(bams)

# A single BAM file
res <- gatk(analysisType="CountReads", pathnameR=pathnameFA,
            pathnameI=pathnameBAMs[1L], verbose=TRUE)
print(res)

# Multiple BAM files
res <- gatk(analysisType="CountReads", pathnameR=pathnameFA,
            pathnameI=pathnameBAMs, verbose=TRUE)
print(res)

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-04-13
# o Created.
############################################################################

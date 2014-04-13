library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "gatk")
if (fullTest) {

bin <- findGATK()
print(bin)

# Setup (writable) local data directory structure
setupExampleData()

dataset <- "GATKResourceBundle";
organism <- "GATKExample";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA and BAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "organisms", organism)
pathnameFA <- Arguments$getReadablePathname("exampleFASTA.fasta", path=path)
print(pathnameFA)

path <- file.path("bamData", dataset, organism)
pathnameBAM <- Arguments$getReadablePathname("exampleBAM.bam", path=path)
print(pathnameBAM)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Call GATK CountReads, cf. '(howto) Run the GATK for the first time'
# http://gatkforums.broadinstitute.org/discussion/1209/howto-run-the-gatk-for-the-first-time
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- gatk(analysisType="CountReads", pathnameR=pathnameFA,
            pathnameI=pathnameBAM, verbose=TRUE)
print(res)

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-04-13
# o Created.
############################################################################

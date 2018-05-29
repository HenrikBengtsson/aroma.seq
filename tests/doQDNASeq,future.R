library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && all(isCapableOf(aroma.seq, c("bwa", "picard")))
fullTest <- fullTest && isPackageInstalled("QDNAseq")
fullTest <- fullTest && isDirectory("annotationData,aroma.seq,private");
fullTest <- fullTest && isDirectory("fastqData,aroma.seq,private");
if (fullTest) {

library("future")
strategies <- c("lazy", "eager")
if (future::supportsMulticore()) strategies <- c(strategies, "multicore")
if (require(pkg <- "future.BatchJobs", character.only=TRUE)) {
  strategies <- c(strategies, "batchjobs_local")
}

dataSet <- "AlbertsonD_2012-SCC,AB042"
organism <- "Homo_sapiens"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism, prefix="human_g1k_v37")
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathR <- "fastqData,aroma.seq,private";
##fqs <- IlluminaFastqDataSet$byName(dataSet, organism=organism, paths=pathR)
path <- IlluminaFastqDataSet$findByName(dataSet, organism=organism, paths=pathR)
fqs <- IlluminaFastqDataSet$byPath(path)
fqs <- extract(fqs, 1:2)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QDNAseq on FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (strategy in strategies) {
  plan(strategy)
  print(plan())

  cns <- doQDNAseq(fqs, reference=fa, binWidth=100, tags=c("*", strategy), verbose=-20)
  print(cns)
}

# Display individual output files
for (ii in seq_along(cns)) print(cns[[ii]])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QDNAseq on BAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- BamDataSet$byName(dataSet, tags="bwa,is,-dups", organism=organism)
print(bams)

# QDNAseq on a single BAM file
bam <- bams[[1]]
print(bam)
cn <- doQDNAseq(bam, binWidth=100, verbose=-20)
print(cn)

for (strategy in strategies) {
  plan(strategy)
  print(plan())

  # QDNAseq on a BAM file set
  cns <- doQDNAseq(bams, binWidth=100, tags=c("*", strategy), verbose=-20)
  print(cns)
}

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-31
# o Created.
############################################################################

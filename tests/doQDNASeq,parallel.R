library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "BiocParallel")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
fullTest <- fullTest && isCapableOf(aroma.seq, "picard")
fullTest <- fullTest && isPackageInstalled("QDNAseq")
fullTest <- fullTest && isDirectory("annotationData,aroma.seq,private");
fullTest <- fullTest && isDirectory("fastqData,aroma.seq,private");
fullTest <- fullTest && isPackageInstalled("BatchJobs")
if (fullTest) {

dataSet <- "AlbertsonD_2012-SCC,AB042"
organism <- "HomoSapiens"

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
cns <- doQDNAseq(fqs, reference=fa, binWidth=100, tags=c("*", "parallel"), verbose=-20)
print(cns)

# Display individual output files
for (ii in seq_along(cns)) print(cns[[ii]])


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QDNAseq on BAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- BamDataSet$byName(dataSet, tags="bwa,is,-dups", organism=organism)
print(bams)

# QDNAseq on a single BAM file
bf <- getFile(bams, 1)
print(bf)
cn <- doQDNAseq(bf, binWidth=100, verbose=-20)
print(cn)

# QDNAseq on a BAM file set
cns <- doQDNAseq(bams, binWidth=100, tags=c("*", "parallel"), verbose=-20)
print(cns)

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-31
# o Created.
############################################################################

library("aroma.seq")
fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
fullTest <- fullTest && isCapableOf(aroma.seq, "picard")
fullTest <- fullTest && isPackageInstalled("qdnaseq")
fullTest <- fullTest && isDirectory("annotationData,aroma.seq,private");
fullTest <- fullTest && isDirectory("fastqData,aroma.seq,private");
if (fullTest) {

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData,aroma.seq,private", "organisms", "Human")
fa <- FastaReferenceFile("human_g1k_v37.fasta", path=path)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData,aroma.seq,private", "AlbertsonD_2012-SCC,AB042", "Generic")
fqs <- FastqDataSet$byPath(path, pattern="[.](fq|fastq)$")
fqs <- extract(fqs, 1:2)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QDNASeq on FASTQ files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cns <- doQDNASeq(fqs, reference=fa, binWidth=100, verbose=-20)
print(cns)

# Display individual BAM files
for (ii in seq_along(cns)) {
  cn <- getFile(cns, ii)
  print(cn)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# QDNASeq on BAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "bamData/AlbertsonD_2012-SCC,AB042,bwa,is,-dups/Generic/"
bams <- BamDataSet$byPath(path)
print(bams)

# QDNASeq on a single BAM file
bf <- getFile(bams, 1)
print(bf)
cn <- doQDNASeq(bf, binWidth=100, verbose=-20)
print(cn)

# QDNASeq on a BAM file set
cns <- doQDNASeq(bams, binWidth=100, verbose=-20)
print(cns)

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-08-31
# o Created.
############################################################################

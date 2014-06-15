library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "YeastTest"
organism <- "SaccharomycesCerevisiae"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# FASTQ data
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bams <- doBowtie2(fqs, reference=fa, verbose=-20)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count nucleotides
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bam <- bams[[1]]
print(bam)

chrs <- getTargetNames(bam)
loci <- data.frame(chromosome=chrs[1L], pos=1:100)
counts <- countNucleotides(bam, loci=loci, verbose=-10)
print(counts)

} # if (fullTest)


############################################################################
# HISTORY:
# 2014-06-15
# o Created.
############################################################################

library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bowtie2")
if (fullTest) {

## Get BAM files
setupExampleData()
dataSet <- "TopHat-example"
organism <- "LambdaPhage"
fa <- FastaReferenceFile$byOrganism(organism)
fqs <- FastqDataSet$byName(dataSet, organism=organism)
bams <- doBowtie2(fqs, reference=fa, verbose=-20)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BAM set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print(bams)

bamsT <- splitByTargetType(bams)
str(bamsT)

bamsT <- splitByTargetType(bams, as="index")
str(bamsT)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BAM file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bam <- bams[[1]]
print(bam)

fs <- getFlagStat(bam)
print(fs)

} # if (fullTest)

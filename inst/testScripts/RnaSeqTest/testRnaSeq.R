# source("testRnaSeq.R")
# Test aroma.seq RNA-seq (through htseq-count)
# (201402 - Need to figure out how to call TopHat2Alignment properly now)

# Cf. system.file("testScripts/RnaSeqProto/rnaSeqScript.R", package="aroma.seq")
# which also needs to be modified to reflect the latest arg types.

options(stringsAsFactors=FALSE)
library(aroma.seq)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check pre-requisites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
stopifnot(isCapableOf(aroma.seq, "bowtie2"))
stopifnot(isCapableOf(aroma.seq, "tophat2"))




# - - - - - - - - - -
# Get annotations
# - - - - - - - - - -

# Set up local annotation dir
Organism <- "Saccharomyces_cerevisiae"
annotsPathLocal <- Arguments$getWritablePath(file.path("annotationData", "organisms", Organism))

# Get gtf file
annotsPath <- system.file("exData/annotationData", package="aroma.seq")
gtfFile <- findFiles(file.path(annotsPath, "organisms", Organism), pattern="[.]gtf[.gz]*$")

# Copy gtf file locally
gtfFileLocal <- file.path(annotsPathLocal, basename(gtfFile))
createLink(link=gtfFileLocal, target=gtfFile, skip=TRUE)
if (regexpr("[.]gz$", gtfFileLocal, ignore.case=TRUE) != -1L) {
  gtfFileLocal <- gunzip(gtfFileLocal, skip=TRUE)
}
gtf <- GtfDataFile(gtfFileLocal)
gtf <- Arguments$getInstanceOf(gtf, "GtfDataFile");

# - - - - - - - - - -
# Get data
# - - - - - - - - - -

# Set up local data dir
DataSet <- "YeastTest"
dataPathLocal <- Arguments$getWritablePath(file.path("fastqData", DataSet, Organism))

# Path to sample reads
dataPath <- system.file("exData/fastqData", package="aroma.seq")

# Copy fastq files locally
fastqFiles <- findFiles(dataPath, pattern="[.]fastq[.gz]*$", recursive=TRUE, firstOnly=FALSE)
sapply(fastqFiles, function(f) {
  createLink(link=file.path(dataPathLocal, basename(f)), target=f)
})


# - - - - - - - - - -
# Build ref index
# - - - - - - - - - -
refFile <- findFiles(path=file.path(annotsPath, "organisms", Organism), pattern="*[.]fa[.gz]*$")
refFileLocal <- file.path(annotsPathLocal, basename(refFile))
createLink(link=refFileLocal, target=refFile)
refFa <- FastaReferenceFile(refFileLocal)
refIndex <- buildBowtie2IndexSet(refFa, verbose=-100)


# - - - - - - - - - -
# Set up input FASTQ dataset, with group names
# - - - - - - - - - -
fqDs <- FastqDataSet$byPath(path=dataPathLocal, pattern="[.]fastq[.gz]*$", paired=TRUE)
fqDs <- setFullNamesTranslator(fqDs, function(names, ...) sub("_", ",", names))


# - - - - - - - - - -
# Call TopHat to align reads
# - - - - - - - - - -
ta <- TopHat2Alignment(dataSet=fqDs, groupBy="name", indexSet=refIndex, transcript=gtf)
bams <- process(ta)
# (201402 => Error: Unknown arguments: prefix)
# res <- getOutputDataSet(ta)


# - - - - - - - - - -
# Count reads
# - - - - - - - - - -

# [ TODO:  Sort TopHat bam files here ]

# Convert TopHat accept_hits.bam to sam
obams <- getPathnames(getOutputDataSet(ta))
osams <- sub(".bam$", ".sam", obams)
for (i in seq_along(obams)) {
  samtoolsView(obams[i], osams[i])
}

for (i in seq_along(osams)) {
  samFile <- osams[i]
  htseqCount(samFile=samFile, gfFile=gtfFile, outFile=sub(".sam$", ".count", samFile))
}


cat("done\n")








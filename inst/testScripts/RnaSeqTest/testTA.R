# source("testTA.R")
# - test TopHat2Alignment()

options(stringsAsFactors=FALSE)
library(aroma.seq)

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
Arguments$getWritablePath(annotsPathLocal)
refFileLocal <- file.path(annotsPathLocal, basename(refFile))
createLink(link=refFileLocal, target=refFile)
refFa <- FastaReferenceFile(refFileLocal)
refIndex <- buildBowtie2IndexSet(refFa, verbose=-100)

# - - - - - - - - - -
# Align reads
# - - - - - - - - - -
fqDs <- FastqDataSet$byPath(path=dataPathLocal, pattern="[.]fastq[.gz]*$", paired=TRUE)
fqDs <- setFullNamesTranslator(fqDs, function(names, ...) sub("_", ",", names))
ta <- TopHat2Alignment(dataSet=fqDs, groupBy="name", indexSet=refIndex, transcripts=gtf)
bams <- process(ta)

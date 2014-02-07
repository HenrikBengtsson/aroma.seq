# source("testRnaSeq.R")
# 201402 - Figure out how to call TopHat2Alignment properly now

# Cf. system.file("testScripts/RnaSeqProto/rnaSeqScript.R", package="aroma.seq")
# which also needs to be modified to reflect the latest arg types.

library(aroma.seq)

# Get sample gtf file
annotsPath <- system.file("exData/annotationData", package="aroma.seq")
Organism <- "SaccharomycesCerevisiae"
gtfFile <- findFiles(file.path(annotsPath, "organisms", Organism), pattern="[.]gtf[.gz]*$")

# Get some sample reads
dataPath <- system.file("exData/fastqData", package="aroma.seq")
DataSet <- "YeastTest"

# Set up local aroma dirs
annotsPathLocal <- Arguments$getWritablePath(file.path("annotationData", "organisms", Organism))
dataPathLocal <- Arguments$getWritablePath(file.path("fastqData", DataSet, Organism))

# Copy gtf file locally
gtfFileLocal <- file.path(annotsPathLocal, basename(gtfFile))
copyFile(gtfFile, gtfFileLocal)
if (regexpr("[.]gz$", gtfFileLocal, ignore.case=TRUE) != -1L) {
  gtfFileLocal <- gunzip(gtfFileLocal)
}
transcripts <- GtfDataFile(gtfFileLocal)
transcripts <- Arguments$getInstanceOf(transcripts, "GtfDataFile");

# Copy fastq files locally
fastqFiles <- findFiles(dataPath, pattern="[.]fastq[.gz]*$", recursive=TRUE, firstOnly=FALSE)
sapply(fastqFiles, function(f) {copyFile(f, file.path(dataPathLocal, basename(f)))})
fastqFilesLocal <- findFiles(dataPathLocal, pattern="[.]fastq[.gz]*$", recursive=TRUE, firstOnly=FALSE)
# - not used


######################################
# Confirm how to call TopHat2Alignment now
######################################

fqDs <- FastqDataSet$byPath(path=dataPathLocal, pattern="[.]fastq[.gz]*$", paired=TRUE)
fqDs <- setFullNamesTranslator(fqDs, function(names, ...) sub("_", ",", names))                          

# Build ref index
refFile <- findFiles(path=file.path(annotsPath, "organisms", Organism), pattern="*[.]fa[.gz]*$")
refFileLocal <- file.path(annotsPathLocal, basename(refFile))
copyFile(refFile, refFileLocal)
refFa <- FastaReferenceFile(refFileLocal)
iSet <- buildBowtie2IndexSet(refFa)

ta <- TopHat2Alignment(dataSet=fqDs, groupBy="name", indexSet=iSet, transcript=transcripts)
bams <- process(ta)
# => Error: Unknown arguments: prefix








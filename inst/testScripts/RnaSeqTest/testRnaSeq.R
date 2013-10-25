# source("testRnaSeq.R")
# Test aroma.seq RNA-seq (through htseq-count)

library(aroma.seq)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set up run
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pathExData <- system.file(file.path("exData"), package="aroma.seq")

# Store run parameters in a list
config <- list()
config$bOverwrite <- TRUE

# Metadata
config$datasetName <- "YeastTest"
config$organism <- "SC"

# Dir for reference fasta and gtf files; these will be copied locally
config$pathRef <- file.path(pathExData, "RNA", "annotationData")

# Dir for input fastq files; these will be copied locally
config$pathData <- file.path(pathExData, "RNA", "fastqData")
config$bPairedEnd <- TRUE

# Not yet used:
# config$qualityEncoding <- "illumina" # c("sanger", "solexa", "illumina"); cf. qrqc:ReadSeqFile

# Set up local paths to annotations and input files
pathLocalAnnots <- file.path("annotationData", "organisms", config$organism)
pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
if (file.exists(config$pathRef)) {
  copyDirectory(from=config$pathRef, to=pathLocalAnnots, overwrite=config$bOverwrite)
}
pathLocalData <- file.path("fastqData", config$datasetName, config$organism)
pathLocalData <- Arguments$getWritablePath(pathLocalData)
if (file.exists(config$pathData)) {
  patternData <- "fastq[.gz]*$"
  dataFiles <- findFiles(path=config$pathData, pattern=patternData, firstOnly=FALSE, recursive=TRUE)
  file.copy(from=dataFiles, to=pathLocalData, overwrite=config$bOverwrite)
}

# Gunzip reference / input data if necessary
sapply(findFiles(path=pathLocalAnnots, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=config$bOverwrite)})
sapply(findFiles(path=pathLocalData, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=config$bOverwrite)})

# Setup reference fasta file
refFas <- FastaReferenceSet$byPath(path=pathLocalAnnots, pattern="[.](fa|fasta)$")
refFasta <- getFile(refFas, 1)  # Presuming there is only one reference fasta file

# Setup input FASTQ set
if (!config$bPairedEnd) {
  fqs <- FastqDataSet$byPath(pathLocalData, pattern="[.](fq|fastq)$")   ## This pattern should already be the default...
} else {
  fqs <- FastqDataSet$byPath(pathLocalData, paired=TRUE);
}

# Check here for bowtie2 capability before proceeding
# (bowtie2 is required by TopHat, i.e. required even if reference index already exists)
stopifnot(isCapableOf(aroma.seq, "bowtie2"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build index, align reads to reference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Build bowtie2 index on reference; [ maybe this should be pre-built for package testing? ]
is <- buildBowtie2IndexSet(refFasta, verbose=-10)  # is = 'index set'

# Align input reads using TopHat
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
process(ta, verbose=TRUE)
taout <- getOutputDataSet(ta)  # For now, this is a GenericDataFileSet


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ TODO:  Sort TopHat bam files here ]

# Gene model
gtfFile <- findFiles(path=pathLocalAnnots, pattern="gtf$")[1]

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


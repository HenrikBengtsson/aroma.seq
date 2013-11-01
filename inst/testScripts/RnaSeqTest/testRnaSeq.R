# source("testRnaSeq.R")
# Test aroma.seq RNA-seq (through htseq-count)

library("aroma.seq")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check pre-requisites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
stopifnot(isCapableOf(aroma.seq, "bowtie2"))
stopifnot(isCapableOf(aroma.seq, "tophat2"))


dataSet <- "YeastTest"
organism <- "SC"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create local unzipped copies of data directories
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The directory where all example data files are
path0 <- system.file(file.path("exData"), package="aroma.seq", mustWork=TRUE)

# Annotation data
path <- file.path("annotationData", "organisms", organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)
sapply(GenericDataFileSet$byPath(path, pattern="[.]gz$"), FUN=gunzip)


# FASTQ data
path <- file.path("fastqData", dataSet, organism)
if (!isDirectory(path)) copyDirectory(from=file.path(path0, path), to=path)
sapply(GenericDataFileSet$byPath(path, pattern="[.]gz$"), FUN=gunzip)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
path <- file.path("annotationData", "organisms", organism)
fas <- FastaReferenceSet$byPath(path)
fa <- getFile(fas, 1)  # Presuming there is only one reference fasta file
print(fa)

# FASTQ data
path <- file.path("fastqData", dataSet, organism)
fqs <- FastqDataSet$byPath(path, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TopHat2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, verbose=TRUE)  # is = 'index set'
print(is)

# Align input reads using TopHat
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
process(ta, verbose=TRUE)

res <- getOutputDataSet(ta)


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


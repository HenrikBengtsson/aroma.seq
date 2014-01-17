# source("protoYeast.R")

library("aroma.seq")
library(edgeR)
source("setupConfig.R")
source("setupDirs.R")
source("checkConfig.R")
options(stringsAsFactors=FALSE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set up run
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

path0 <- system.file(file.path("exData"), package="aroma.seq", mustWork=TRUE)
Organism <- "SaccharomycesCerevisiae"
DataSetName <- "YeastTest"
Paths <- list(ref=file.path(path0, "annotationData", "organisms", organism),
              data=file.path(path0, "fastqData", DataSetName, organism))

config <- setupConfig(pathRef=Paths$ref,
                      pathData=Paths$data,
                      datasetName=DataSetName,
                      organism=Organism,
                      bPairedEnd=TRUE,   ## This is actually not used
                      bOverwrite=FALSE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set up dirs under aroma working dir
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
system.time({
  config <- setupDirs(config=config, bGunzip=TRUE)
})
save(config, file="config_000.RData")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set up reference index
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

refFas <- FastaReferenceSet$byPath(path=config$pathLocalAnnots, pattern="[.](fa|fasta)$")
refFasta <- getFile(refFas, 1)  # Presuming there is only one reference fasta file
stIndex <- 
  system.time({
    iSet <- buildBowtie2IndexSet(refFasta, verbose=10)
  })
print(iSet)
save(iSet, file="iSet_000.RData")
save(stIndex, file="stIndex_000.RData")  # i.e. save timing

# Copy reference fasta file to bowtie index directory
# - This is for TopHat; otherwise it builds the reference .fa from the bowtie2 index
file.copy(from=getPathname(refFasta), to=getPath(iSet), overwrite=FALSE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sets up gene model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gtfs <- GenericDataFileSet$byPath(config$pathLocalAnnots, pattern="[.]gtf$")
gtf <- gtfs[[1]]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup input FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Create FastqDataSet
fqsPE <- FastqDataSet$byPath(path=config$pathLocalData, paired=TRUE)
fqsPE <- setFullNamesTranslator(fqsPE, function(names, ...) sub("_", ",", names))
print(getFullNames(fqsPE))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Align reads w/ TopHat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (length(fqsPE) > 0) {
  stDoTopHatPE <- 
    system.time({
      ta <- TopHat2Alignment(fqsPE, groupBy="name", indexSet=iSet, transcripts=gtf)
      bamsPE <- process(ta)
    })
  save(stDoTopHatPE, file="stDoTopHatPE_000.RData")
  save(bamsPE, file="bamsPE_000.RData")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads with htseq
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (exists("bamsSE")) {
  stHCSE <- system.time({
    res <-
      lapply(bamsSE, function(bam, gtfFile) {
        pathname <- getPathnames(bam)
        htseqCount(pathnameS=pathname, gff=gtfFile)
      }, gtfFile=getPathname(gtf)) 
  })
  save(stHCSE, file="stHCSE_000.RData")
}
if (exists("bamsPE")) {
  stHCPE <- system.time({
    res <-
      lapply(bamsPE, function(bam, gtfFile) {
        pathname <- getPathname(bam)
        htseqCount(pathnameS=pathname, gff=gtfFile)
      }, gtfFile=getPathname(gtf)) 
  })
  save(stHCPE, file="stHCPE_000.RData")
}
htCountFiles <- GenericDataFileSet$byPath(path="tophat2Data", pattern="[.]count$", recursive=TRUE)

# Create human-readable count matrix
featureNames <- read.table(getPathnames(htCountFiles)[1], sep="\t")[,1]
countMat <- NULL
for (i in seq_along(htCountFiles)) {
  cat(i, " ")
  df <- read.table(getPathnames(htCountFiles)[i], sep="\t")
  stopifnot(identical(featureNames, df[,1]))
  countMat <- cbind(countMat, df[,2])
}
rownames(countMat) <- featureNames
colnames(countMat) <- sapply(getPathnames(htCountFiles), function(path) {basename(getParent(path))})
countMat <- countMat[,order(colnames(countMat))]
save(countMat, file="countMat_000.RData")
countMatOut <- countMat
countMatOut <- cbind(featureNames, countMat)
colnames(countMatOut) <- c("Feature", colnames(countMat))
write.table(countMatOut, file="countMat_000.xls", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# For edgeR
dge <- readDGE(getPathnames(htCountFiles), columns=c(1,2))
save(dge, file="dge_000.RData")

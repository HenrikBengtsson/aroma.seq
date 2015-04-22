library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {


# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "TopHat-example"
organism <- "LambdaPhage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, verbose=-10)
print(is)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)


bams <- list()

## Destination directory
pathD <- NULL

## For each FASTQ file, split into four, align and merge
size <- 1/4

for (ii in seq_along(fqs)) {
  fq1 <- fqs[[ii]]
  name <- getFullName(fq1)
  
  fq2 <- getMateFile(fq1)
  
  path <- file.path("fastqData", sprintf("tmp-%s", name), organism)

  fqs1M <- splitUp(fq1, size=size, path=path, verbose=TRUE)
  fqs2M <- splitUp(fq2, size=size, path=path, verbose=TRUE)
  fqsM <- newInstance(fqs, as.list(fqs1M), paired=TRUE)

  # BWA with BWA 'aln' options '-n 2' and '-q 40'.
  alg <- BwaAlignment(fqsM, indexSet=is, n=2, q=40, verbose=-10)
  print(alg)

  ## Infer output directory
  if (is.null(pathD)) {
    tags <- c(getTags(alg), "split")
    datasetD <- paste(c(getFullName(fqs), tags), collapse=",")
    pathD <- file.path(getRootPath(alg), datasetD, getOrganism(alg))
    pathD <- Arguments$getWritablePath(pathD)
  }

  bamsM <- process(alg, verbose=-10)
  print(bamsM)

  ## Merge aligned BAM files
  groupBy <- list(seq_along(bamsM))
  names(groupBy) <- name
  bm <- BamMerger(bamsM, groupBy=groupBy)
  print(bm)
  bamM <- process(bm, verbose=-10)[[1]]
  print(bamM)

  pathM <- getPath(bamM)
  bam <- renameTo(bamM, path=pathD)
  bai <- buildIndex(bam)
  print(bam)

  ## Cleanup
  removeDirectory(getPath(fqs1M), recursive=TRUE)
  removeDirectory(getPath(bamsM), recursive=TRUE)
  removeDirectory(pathM, recursive=TRUE)
  
  bams[[ii]] <- bam
} # for (ii ...)

bams <- BamDataSet(bams)
print(bams)


} # if (fullTest)



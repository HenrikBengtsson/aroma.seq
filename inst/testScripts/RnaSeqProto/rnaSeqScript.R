# 
# Prototype RNA-seq script
# 
# This script is a pipeline for performing RNA-seq read quantification in aroma.seq.
# It starts with fastq files and ends with a read count matrix on transcript features.
# Sample code to set up the aroma.seq package and other suggested R packages, and download
# reference genome and transcript annotations are included for convenience.
# External binaries bowtie2 and tophat must exist on the user's path (setup described elsewhere).
# 
# This script builds (or assumes) a directory structure as follows:
#   ./annotationData/organisms/<organism>/
#     - This contains (or will contain) reference annotation .fa, .gtf, bowtie2/bwa indices, etc.
#   ./fastqData/<dataset>/<organism>/, which contains input *.fastq files
#     - This contains the reads that will be aligned/quantified.  Currently all sample reads must be paired-end or single-end, not mixed.
#
# Output dir at end of pipeline:  tophat2Data/ (htseq-count output goes here as well)
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup R and run config
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Set up R packages
bOverwrite <- TRUE
bSetupR <- FALSE
if (bSetupR) {
  # Setup BioConductor / aroma packages
  source( "http://www.bioconductor.org/biocLite.R" )
  biocLite("BiocUpgrade")
  biocLite( c("ShortRead","DESeq", "edgeR") )
  source("http://www.braju.com/R/hbLite.R")
  hbInstall("aroma.seq", devel=TRUE)
}

# Set up run parameters
config <- list()

# Dir for reference fasta and gtf files; these will be copied locally
config$pathRef <- "" # E.g. "/data/annotationData/organisms/HomoSapiens"

# Dir for input fastq files; these will be copied locally
config$pathData <- "" # E.g. "/data/SRA/GSE18508/"

# Dataset metadata
config$datasetName <- "MyDataSet"
config$organism <- "HomoSapiens"
config$bPairedEnd <- FALSE
config$qualityEncoding <- "illumina" # c("sanger", "solexa", "illumina"); cf. qrqc:ReadSeqFile
# - This is difficult to determine automatically or for users to obtain, but important
#   to get right for the quality assessment step.  In particular, the code bombs if it
#   encounters phred scores outside range (though it might be worse if it does not bomb).


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Start aroma.seq run
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library("aroma.seq")
pathExData <- system.file("exData", package="aroma.seq")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Copy reference and input files to local dirs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# If config$pathRef and config$pathData do not exist, then the reference annotations
# and input fastq files must already be in the correct local locations.
pathLocalAnnots <- file.path("annotationData", "organisms", config$organism)
pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
if (exists("config$pathRef")) {
  copyDirectory(from=config$pathRef, to=pathLocalAnnots, overwrite=bOverwrite)
}
pathLocalData <- file.path("fastqData", config$datasetName, config$organism)
pathLocalData <- Arguments$getWritablePath(pathLocalData)
if (exists("config$pathData")) {
  patternData <- "fastq$"
  dataFiles <- findFiles(pattern=patternData, path=pathData, firstOnly=FALSE)
  file.copy(from=dataFiles, to=pathLocalData, overwrite=bOverwrite)
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Gunzip input / reference data if necessary
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sapply(findFiles(path=pathLocalData, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=bOverwrite)})
sapply(findFiles(path=pathLocalAnnots, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=bOverwrite)})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup reference fasta file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

refFas <- FastaReferenceSet$byPath(path=pathLocalAnnots, pattern="[.](fa|fasta)$")
# Download reference files, gunzip, move to standard location, if needed
if (length(refFas) < 1)
{
  
  # Get list of reference URLs from the package
  source(file.path(pathExData, "ReferenceGenomes.R"))
  
  RefUrls <- RefUrlsList[[config$organism]]
  sapply(RefUrls, function(loc)
  {
    fnameGZ <- basename(loc)
    fname <- sub(".gz", "", fnameGZ)   # [ Cludgey on two counts: assuming .gz suffix, assuming gzip'd in the first place ]
    # TODO:  Make this agnostic to .gz or .bz2 or uncompressed
    
    # Download file if it has not been downloaded already
    if (!isFile(file.path(pathLocalAnnots, fnameGZ)) && !isFile(file.path(pathLocalAnnots, fname)))
    {
      pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
      downloadFile(loc, path=".")
      # gunzip(fnameGZ)  ## most indexers can handle .gz input, so skip this
      renameFile(fnameGZ, file.path(pathLocalAnnots, fnameGZ))
      returnVal <- file.path(pathLocalAnnots, fnameGZ)
    } else {
      returnVal <- NULL
    }
  })
}
## [ LIMITATION: Presuming there is only one reference fasta file ]
refFasta <- getFile(refFas, 1)
print(refFasta)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup input FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (!bPairedEnd) {
  fqs <- FastqDataSet$byPath(pathLocalData, pattern="[.](fq|fastq)$")   ## This pattern should already be the default...
  print(fqs)
} else {
  fqs <- FastqDataSet$byPath(pathLocalData, paired=TRUE);
  print(fqs)
}

# Check here for bowtie2 capability before proceeding
# (bowtie2 is required by TopHat, i.e. required even if reference index already exists)
stopifnot(isCapableOf(aroma.seq, "bowtie2"))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set up reference index
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

is <- buildBowtie2IndexSet(refFasta, verbose=-10)  # is = 'index set'
print(is)
# - NB bowtie2 is actually only run of the index does not exist already


# Copy reference fasta file to bowtie index directory
# (This is for TopHat; otherwise it builds the reference .fa from the bowtie2 index)
copyFile(getPathname(refFasta), file.path(dirname(getIndexPrefix(is)), getFilename(refFasta)),
         overwrite=bOverwrite)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Quality assess
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (FALSE) {  # Skip for now
  if (!bPairedEnd)
  {
    # NB: For now, the following overwrites the fastq dataset
    
    # adapter trim first
    fqs <- doScythe(fqs, "pairedEnd")
    
    # quality trim reads
    fqs <- doSickle(fqs, "pairedEnd")
    
  } else {
    pdfList <- sapply(fqs, report)    ### Uses qrqc; this could take a while...probably should make verbose
    fqs <- sapply(fqs, doScythe)
    fqs <- sapply(fqs, doSickle)
  }
  # Quality assessment reports
  pdfs <- report(fqs)
  # doQRQC(fqs)   ## TBD
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Align reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Gene model
gtfFile <- findFiles(path=pathLocalAnnots, pattern="gtf$")[1]
# - POSSIBLE MANUAL STEP:  User needs to specify which gene model to use if more than one is available

ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
# TODO:  Make this run only if not done already
system.time({
  process(ta, verbose=TRUE)
})
taout <- getOutputDataSet(ta)  # For now, this is a GenericDataFileSet

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Organize, sort and index the BAM files and create SAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ ... not done yet ... ]
# [ This step needed for htseq-count; use Rsamtools as alternative? ]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (FALSE) {
  # Convert bam to sam
  ods <- getOutputDataSet(ta)
  acceptedHits <- sapply(sapply(ods, getPathnames, simplify=FALSE), function(x) {x[grep("accepted", x)]})
  acceptedHitsSam <- sapply(acceptedHits, function(f) {Arguments$getWritablePathname(sub(".bam$", ".sam", f))})
  sapply(seq_along(acceptedHits), function(i) {samtoolsView(acceptedHits[i], acceptedHitsSam[i])})
  for (i in seq_along(acceptedHits)) {
    samtoolsView(acceptedHits[i], acceptedHitsSam[i])
  }
  # - ADD TO WISH LIST: Replace acceptedHitsSam w/ in-memory object (say?); i.e. temporary / on-demand rather than
  #   permanent addition to TopHatAlignment / file system.
}

# Convert TopHat accept_hits.bam to sam
obams <- getPathnames(getOutputDataSet(ta))
osams <- sub(".bam$", ".sam", obams)
for (i in seq_along(obams)) {
  samtoolsView(obams[i], osams[i])
}

system.time({
  for (i in seq_along(osams)) {
    samFile <- osams[i]
    htseqCount(samFile=samFile, gfFile=gtfFile, outFile=sub(".sam$", ".count", samFile))
  }
})

# NB: The following is not necessary, since readDGE() does this for us
# Create the count matrix 'by hand', for human purposes.
if (FALSE) {
  # Load the htseq count files as a TabularTextFileSet
  db <- TabularTextFileSet$byPath(path=getPath(ta), pattern="[.]count$", recursive=TRUE)
  ## [ ISSUE:  This drops the first row, since treated as column names ]
  
  # Get the 2nd column from each htseq count table
  mat <-
    sapply(seq_along(db), function(i) {
      ttf <- getFile(db, i)
      col2 <- readColumns(ttf,2)
    }, simplify=TRUE)
  mat <- sapply(mat, function(x) x)  # Silly but works to convert list to matrix (...)
  rns <- readColumns(getFile(db, 1), 1)[,1]
  rownames(mat) <- rns
  bDrop <- rns %in% c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique")
  mat <- mat[!bDrop,]
}


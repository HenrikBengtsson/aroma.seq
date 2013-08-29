# 201307-08 Taku
# Prototype RNA-seq workflow (minimal checking)
#
# This script builds / assumes a directory structure as follows (20130730 convention):
#   ./annotationData/organisms/<organism>/<reference annotation .fa, .gtf, bowtie2 / bwa indices, etc.>
#   ./fastqData/<dataset>/<organism>/*.fastq
# 
# Output dirs:  probably bamData/  (e.g. output from BWA, bowtie2)
# - TopHat output probably goes into its own output folder
# 

bDEBUG <- FALSE

library("aroma.seq")
## [- NB:  Order of library loading may be significant for the code to work properly ]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup config variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ Using full paths here only for local 'dump' dirs ]
# If pathRef does not exist, then reference genome will be downloaded
pathRef <- "/Users/tokuyasu/Data/CCSP/annotationData/organism/SaccharomycesCerevisiae"
pathData <- "/Users/tokuyasu/Data/CCSP/20130620.Yeast"
datasetName <- "YeastTest"
organism <- "SaccharomycesCerevisiae"
bPairedEnd <- TRUE
overwrite <- FALSE
qualityEncoding <- "illumina" # c("sanger", "solexa", "illumina"); cf. qrqc:ReadSeqFile
# - This is difficult to determine automatically, and difficult for naive users to obtain, but important
#   to get right for the quality assessment step.  In particular, the code bombs if it encounters phred
#   scores outside range (might be worse if it does not bomb actually).
bSetupR <- FALSE

# For edgeR
MinNumberOfReplicates <- 2  # Used in edgeR preprocessing as recommended by Anders et al; NB must be >= 3 for edgeR to work, I think
Groups <- c("a", "a", "b", "b")  # Sample assignments to (two) groups


# R packages, reference genome downloads
if (bSetupR) {
  # Setup BioConductor packages
  source( "http://www.bioconductor.org/biocLite.R" )
  biocLite("BiocUpgrade")
  biocLite( c("ShortRead","DESeq", "edgeR") )
  
  ## [ Henrik:  IS THIS THE BEST WAY TO SETUP AROMA.SEQ? ]
  source("http://aroma-project.org/hbLite.R")
  hbLite("aroma.seq")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup local directories for annotations and data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ Copying these to a local dir is a temp step; needed for now to handle read-only dirs
#   that don't accept output (aroma convention is to create local dirs for the output only?) ]
pathLocalAnnots <- file.path("annotationData", "organisms", organism)
pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
if (exists("pathRef"))
{
  copyDirectory(from=pathRef, to=pathLocalAnnots, overwrite=FALSE)  
}
pathLocalData <- file.path("fastqData", datasetName, organism)
pathLocalData <- Arguments$getWritablePath(pathLocalData)
# patternData <- "fastq$"
patternData <- "*1e\\+.*fastq$"  # subsampled data
dataFiles <- findFiles(pattern=patternData, path=pathData, firstOnly=FALSE)
file.copy(from=dataFiles, to=pathLocalData, overwrite=FALSE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ LIMITATION:  This will fail if both the .gz and the gunzip'd reference file exist in pathLocalAnnots ]
refFas <- FastaReferenceSet$byPath(path=pathLocalAnnots, pattern="[.](fa|fasta)(.gz)*$")
# Download reference fasta file if necessary
if (length(refFas) < 1)
{
  # Download reference files, gunzip, move to standard location, if needed
  load("RefUrlsList.RData")
  RefUrls <- RefUrlsList[[organism]]
  sapply(RefUrls, function(loc)
  {
    fnameGZ <- basename(loc)
    fname <- sub(".gz", "", fnameGZ)   # [ Cludgey on two counts: assuming .gz suffix, assuming gzip'd in the first place ] 
    # TODO:  Make this agnostic to .gz or .bz2 or not
    
    # Check to see if the file has been downloaded already
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
  if (FALSE)  {# First version:  paired end data as a list
    fqs1 <- FastqDataSet$byPath(pathLocalData, pattern="[._]1.(fq|fastq)$")  ## (Should allow these to be .gz as wel)
    fqs2 <- FastqDataSet$byPath(pathLocalData, pattern="[._]2.(fq|fastq)$")
    fqs <- list(read1 = fqs1, read2 = fqs2)   ## [ Is this reasonable?  fqs = list => PE data? ]
  } else { # New paired-end style FastqDataSet
    fqs <- FastqDataSet$byPath(pathLocalData, paired=TRUE);
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 0: Set up reference index
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ LIMITATION:  Following should not be req'd if reference indices already exist ]
isCapableOf(aroma.seq, "bowtie2")


# NB: bowtie2 is run only if the index set does not already exist
system.time({
  is <- buildBowtie2IndexSet(refFasta, verbose=-10)  # is = 'index set'
  print(is)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 1: Quality assess
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
# Step 2: Gather experiment metadata
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [TBD:  This is probably a manual step for the user;
#  important especially if the input filenames are not sufficient / satisfactory for the user ]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 3: Align reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Gene model
gtfFile <- findFiles(path=pathLocalAnnots, pattern="gtf$", firstOnly=TRUE)
# - MANUAL STEP:  User must specify which gene model to use, if more than one is available

if (FALSE)  # Example of how to run tophat 'by hand'
{
  fqs1 <- extract(fqs, 1)
  Reads1 <- getPathname(getFile(fqs1, 1))
  Reads2 <- getPathname(getMateFile(getFile(fqs1,1)))
  # Alt1:
  # Reads1 <- sapply(getFilePairs(fqs1)[,1], function(x) {getPathname(x)})[1]
  # Reads2 <- sapply(getFilePairs(fqs1)[,2], function(x) {getPathname(x)})[1]
  # Alt2:
  # Reads1 <- getPathnames(fqs1)[1]
  # Reads2 <- getPathname(getMateFile(getFile(fqs1, 1)))    
  # source("/Users/tokuyasu/SVN/aroma.seq.RS/R/tophat.R")
  # source("/Users/tokuyasu/SVN/aroma.seq.RS/R/systemTopHat.R")
  ## - Possible remaining bug for bin assignment => use Sys.which("tophat")
  system.time({
    tophat(getIndexPrefix(is), reads1=Reads1, reads2=Reads2, verbose=TRUE)
  })
  system.time({
    tophat(getIndexPrefix(is), reads1=Reads1, reads2=Reads2, optionsVec=c("G"=gtfFile, "o"="tophat_byHand"), verbose=TRUE)
  })
}

source("/Users/tokuyasu/work/projects/12/CCSP/proto/TopHat2Alignment.R")
if (FALSE) {  # This is DEBUG-only
  source("/Users/tokuyasu/work/projects/12/CCSP/proto/tophat.R")
  source("/Users/tokuyasu/work/projects/12/CCSP/proto/systemTopHat.R")  ## This is key - need to be able to find the tophat executable!
}
ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
# Following is DEBUG-only
if (bDEBUG) {
  ta <- TopHat2Alignment(dataSet=extract(fqs, files=c(1,2)), indexSet=is)
}

# TODO:  Make this run only if not done already
system.time({
  process(ta, verbose=TRUE)
})

taout <- getOutputDataSet(ta)  # For now, this is a GenericDataFileSet

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 4: Organize, sort and index the BAM files and create SAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ ... I'm not sure this step is really necessary; probably is needed for htseq-count; use an Rsamtools alternative? ] ;

## [ ...Step 5 is now Step 8... ]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 6: Count reads
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

# Creating the count matrix 'by hand' may be useful, but skip for now, since readDGE() does this for us
if (FALSE) {
  # Load the htseq count files as a TabularTextFileSet
  db <- TabularTextFileSet$byPath(path=getPath(ta), pattern="[.]count$", recursive=TRUE)
  # - ISSUE:  This drops the first row, since treated as column names
  
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 7a:  Set up edgeR count datastruct and filter reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Step 1 edgeR Create container for count data and filter features ;
## Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count:

library(edgeR);
hcounts <- sub(".sam$", ".count", osams)
counts <- readDGE(hcounts)$counts;

## In edgeR, we recommend removing features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates (here, n=3 for the Knockdown group). Filter these as well as non-informative (e. g., non-aligned) features using a command like: ;
noint <- rownames(counts) %in%
  c("no_feature","ambiguous","too_low_aQual",
    "not_aligned","alignment_not_unique")
cpms <- cpm(counts);

# From Anders et al
## 'In edgeR, we recommend removing features without at least 1 read per million in n of
##  the samples, where n is the size of the smallest group of replicates'
keep <- (rowSums(cpms>1)>=MinNumberOfReplicates) & !noint;
counts <- counts[keep,];



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 7b: Visualize and create DGEList for edgeR
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Visualize and inspect the count table using: ;


colnames(counts) <- sub("_1$", "", getFullNames(getInputDataSet(ta)))
# head( counts[,order(samples$condition)], 5 );

## Create a DGEList object (s container for RNA-seq count data), as follows: ;
d <- DGEList(counts=counts, group=Groups);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 8 (was Step 5): Inspect alignments with IGV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ ... This is probably a useful but manual step; user intervention required; skip for auto-pipeline ] ;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run edgeR for differential expression
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# [...See Adam...] ;
#



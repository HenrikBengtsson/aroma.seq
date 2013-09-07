# 201307-09 Taku
# Prototype RNA-seq workflow (minimal checking)
#
# This script builds (or assumes) a directory structure as follows (20130730 convention):
#   ./annotationData/organisms/<organism>/<reference annotation .fa, .gtf, bowtie2 / bwa indices, etc.>
#   ./fastqData/<dataset>/<organism>/*.fastq
# 
# Output dir:  tophat2Data/ (htseq-count output goes here as well)
# 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup config variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

bTest <- TRUE  ## Run script in package test mode
bDebug <- TRUE
bSetupR <- FALSE
bOverwrite <- TRUE

if (bSetupR) {
  # Setup BioConductor / aroma packages
  source( "http://www.bioconductor.org/biocLite.R" )
  biocLite("BiocUpgrade")
  biocLite( c("ShortRead","DESeq", "edgeR") )
  source("http://www.braju.com/R/hbLite.R")
  hbInstall("aroma.seq")
}

library("aroma.seq")

if (bTest) {
  pathD <- system.file("exData", package="aroma.seq")
  for (dir in c("annotationData", "fastqData")) {
    copyDirectory(file.path(pathD, "RNAseq", dir), to=dir, overwrite=bOverwrite)
  }
} else {
  # NB:  Using full paths below.  The following is meant as prototype code to aid a user who has a set
  # of fastq files in pathData, and possibly some reference annotation data in pathRef as well, and
  # wants to get started.  If pathRef does not exist, then reference genome for 'organism' will be
  # downloaded.
  pathRef <- "/Users/tokuyasu/Data/CCSP/annotationData/organisms/SaccharomycesCerevisiae"
  pathData <- "/Users/tokuyasu/Data/CCSP/20130620.Yeast"
}
organism <- "SaccharomycesCerevisiae"
datasetName <- "YeastTest"
bPairedEnd <- TRUE
qualityEncoding <- "illumina" # c("sanger", "solexa", "illumina"); cf. qrqc:ReadSeqFile
# - This is difficult to determine automatically, and difficult for naive users to obtain, but important
#   to get right for the quality assessment step.  In particular, the code bombs if it encounters phred
#   scores outside range (might be worse if it does not bomb actually).

# Sample config for edgeR (not run here)
MinOKReplicates <- 2  # Min number of OK replicate samples per feature; recommended by Anders et al
                      # NB: edgeR requires at least 3 replicate samples to work, I think
Groups <- c("a", "a", "b", "b")  # Sample assignments to (two) groups


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup local directories for annotations and data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ Copying these to a local dir is a temp step; needed for now to handle read-only dirs
#   that don't accept output (aroma convention is to create local dirs for the output only?) ]
pathLocalAnnots <- file.path("annotationData", "organisms", organism)
pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
if (exists("pathRef")) {
  copyDirectory(from=pathRef, to=pathLocalAnnots, overwrite=bOverwrite)  
}
if (exists("pathData")) {
  pathLocalData <- file.path("fastqData", datasetName, organism)
  pathLocalData <- Arguments$getWritablePath(pathLocalData)
  # patternData <- "fastq$"
  patternData <- "*1e\\+.*fastq$"  # subsampled data
  dataFiles <- findFiles(pattern=patternData, path=pathData, firstOnly=FALSE)
  file.copy(from=dataFiles, to=pathLocalData, overwrite=bOverwrite)
} else {
  pathLocalData <- dirname(findFiles(path="fastqData", pattern="[.](fq|fastq)(.gz)*$", recursive=TRUE))
}

# Gunzip input / reference data
sapply(findFiles(path=pathLocalData, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=bOverwrite)})
sapply(findFiles(path=pathLocalAnnots, pattern=".gz$", firstOnly=FALSE),
       function (f) {gunzip(f, overwrite=bOverwrite)})

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

# [ LIMITATION:  Following should not be req'd if reference index already exists ]
isCapableOf(aroma.seq, "bowtie2")


# NB: bowtie2 is run only if the index set does not already exist
system.time({
  is <- buildBowtie2IndexSet(refFasta, verbose=-10)  # is = 'index set'
  print(is)
})

# For TopHat convenience (otherwise it builds reference .fa from bowtie2 index)
copyFile(getPathname(refFasta), file.path(dirname(getIndexPrefix(is)), getFilename(refFasta)),
         overwrite=bOverwrite)


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
gtfFile <- findFiles(path=pathLocalAnnots, pattern="gtf$")[1]
# - POSSIBLE MANUAL STEP:  User needs to specify which gene model to use if more than one is available

ta <- TopHat2Alignment(dataSet=fqs, indexSet=is)
# TODO:  Make this run only if not done already
system.time({
  process(ta, verbose=TRUE)
})
taout <- getOutputDataSet(ta)  # For now, this is a GenericDataFileSet

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 4: Organize, sort and index the BAM files and create SAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# [ This step may be needed for htseq-count, esp. for non-TopHat alignment; use Rsamtools as alternative? ]


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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 7a:  Set up edgeR count datastruct and filter reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Step 1 edgeR Create container for count data and filter features ;
## Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count:

library(edgeR);
hcounts <- sub(".sam$", ".count", osams)
counts <- readDGE(hcounts)$counts;

# Filter features w/ too low coverage and uninformative read counts 
noint <- rownames(counts) %in%
  c("no_feature","ambiguous","too_low_aQual",
    "not_aligned","alignment_not_unique")
cpms <- cpm(counts);
# From Anders et al: 'In edgeR, we recommend removing features without at least 1 read per million in n of
# the samples, where n is the size of the smallest group of replicates'
if (!bTest) {
  keep <- (rowSums(cpms>1)>=MinOKReplicates) & !noint;
  counts <- counts[keep,];
}

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
#
# [ Example code ]
# # Estimate normalization factors using:
# d <- calcNormFactors(d)
# 
# # For simple designs, estimate tagwise dispersion estimates using:
# d <- estimateCommonDisp(d)
# d <- estimateTagwiseDisp(d)
# 
# # For a simple two-group design, perform an exact test for the difference in expression
# # between the two conditions:
# de <- exactTest(d, pair=c("CTL","KD"))
# 



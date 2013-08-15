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
  fqs1 <- FastqDataSet$byPath(pathLocalData, pattern="[._]1.(fq|fastq)$")  ## (Should allow these to be .gz as wel)
  fqs2 <- FastqDataSet$byPath(pathLocalData, pattern="[._]2.(fq|fastq)$")
  print(fqs1)
  print(fqs2)
  fqs <- list(read1 = fqs1, read2 = fqs2)   ## [ Is this reasonable?  fqs = list => PE data? ]
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 0: Set up reference index
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ LIMITATION:  Following should not be req'd if reference indices already exist ]
isCapableOf(aroma.seq, "bowtie2")

## TODO:  Run only if index set does not exist
system.time({
  is <- buildBowtie2IndexSet(refFasta, verbose=-10)  # is = 'index set'
  print(is)
})
# user  system elapsed 
# 15.998   0.310  16.946 


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 1: Quality assess
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
# pdfs <- report(fqs)
# doQRQC(fqs)   ## TBD


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 2: Gather experiment metadata
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# ... [TBD; probably a manual step; how to capture the experimental design?]


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 3: Align reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

source("/Users/tokuyasu/SVN/aroma.seq.RS/R/AbstractAlignment.R")
# - This is Taku's version, which allows dataset to be a list of FastqDataSet's


if (FALSE)  # run tophat 'by hand' first
{
  Reads1 <- getPathnames(fqs[[1]])[1]
  Reads2 <- getPathnames(fqs[[2]])[1]
  source("/Users/tokuyasu/SVN/aroma.seq.RS/R/tophat.R")
  source("/Users/tokuyasu/SVN/aroma.seq.RS/R/systemTopHat.R")
  ## - Crap, the needed fix is to make bin non-NULL using Sys.which("tophat")
  system.time({
    tophat(getIndexPrefix(is), reads1=Reads1, reads2=Reads2, verbose=TRUE)
  })
}


ta <- TopHatAlignment(dataSet=fqs, indexSet=is)
# - Works! (with Taku's version of AbstractAlignment)
process(ta)


stop("Got to here!")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 4: Organize, sort and index the BAM files and create SAM files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ ... I'm not sure this step is really necessary; probably is needed for htseq-count; use an Rsamtools alternative? ] ;

## [ ...Step 5 is now Step 8... ]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 6: Count reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# proposed code ;
counts <- countReads(pBAM(s), geneModel=..., method='htseq');
# - ( Presumably 'counts' should be a class object as well? At least store location of any file output. )

## HtseqDataSet, HtseqDataFile
## - inherited from Generic...
##  (Does htseq have a suggested way to import in R its export data?)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 7a:  Set up edgeR count datastruct and filter reads
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Step 1 edgeR Create container for count data and filter features ;
## Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count:

library("edgeR");
counts <- readDGE(samples$countf)$counts;

## In edgeR, we recommend removing features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates (here, n=3 for the Knockdown group). Filter these as well as non-informative (e. g., non-aligned) features using a command like: ;
noint <- rownames(counts) %in%
  c("no_feature","ambiguous","too_low_aQual",
    "not_aligned","alignment_not_unique")
cpms <- cpm(counts);

# [ The following command might be messed up ]
keep <- (rowSums(cpms>1)>=3) && !noint;
counts <- counts[keep,];

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 7b: Visualize and create DGEList for edgeR
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Visualize and inspect the count table using: ;
colnames(counts) <- samples$shortname;
head( counts[,order(samples$condition)], 5 );

## Create a DGEList object (s container for RNA-seq count data), as follows: ;
d <- DGEList(counts=counts, group=samples$condition);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Step 8 (was Step 5): Inspect alignments with IGV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# [ ... This is probably a useful but manual step; user intervention required; skip for auto-pipeline ] ;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run edgeR for differential expression
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# [...See Adam...] ;
#



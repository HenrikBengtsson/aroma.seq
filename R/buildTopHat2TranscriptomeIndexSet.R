###########################################################################/**
# @set "class=Bowtie2IndexSet"
# @RdocMethod buildTopHat2TranscriptomeIndexSet
# @alias buildTopHat2TranscriptomeIndexSet
#
# @title "Calls TopHat to build a transcriptome index; 'this' is the reference genome index set"
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{gtf}{GtfDataFile to be indexed}
#   \item{outPath}{(optional) Output directory for index and log file.}
#   \item{tiPrefix}{(optional) Prefix for transcriptome index.}
#   \item{...}{Arguments passed to tophat().}
#   \item{skip}{If @TRUE, the index files are not rebuilt if already available.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("buildTopHat2TranscriptomeIndexSet", "Bowtie2IndexSet", function(this,
                                                                             gtf,
                                                                             outPath=NULL,
                                                                             tiPrefix=NULL,
                                                                             ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument gtf
  gtfFile <- getPathname(gtf)
  gtfFile <- Arguments$getReadablePathname(gtfFile)
  
  # Argument outPath
  if (is.null(outPath)) {
    outPath <- file.path(getPath(gtf), "tophat2", getFullName(gtf))
  }
  Arguments$getWritablePath(outPath)  
  
  # Argument tiPrefix
  if (is.null(tiPrefix)) {
    tiPrefix <- "."
  }
  
  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose), add=TRUE);
  }

  # Check for tophat version
  verT <- attr(findTopHat2(), "version");
  if (verT < '2.0.10') {
    throw("TopHat version >= 2.0.10 required")
  }

  # (Pre-existing) index for the reference genome
  refIdxPrefix <- getIndexPrefix(this)

  # Check for existing transcriptome index files
  res <- tryCatch({
    Bowtie2IndexSet$byPrefix(file.path(outPath, tiPrefix, tiPrefix))
  }, error=function(ex) Bowtie2IndexSet());

  # Nothing todo?
  if (skip && isComplete(res)) {
    verbose && cat(verbose, "Transcriptome indexing already done. Skipping.");
    return(res)
  }

  # Call TopHat executable
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  verbose && enter(verbose, "Building transcriptome index.");
  res <- tophat(refIdxPrefix, gtf=gtfFile, outPath=outPath, optionsVec=c("--transcriptome-index"=tiPrefix), ...)
  
  # Locate index set to return
  res <- tryCatch({
    Bowtie2IndexSet$byPrefix(file.path(outPath, tiPrefix, tiPrefix));
  }, error=function(ex) Bowtie2IndexSet());
  
  verbose && exit(verbose);
  
  res
}) # buildTopHat2TranscriptomeIndexSet()

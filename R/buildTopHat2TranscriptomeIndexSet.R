###########################################################################/**
# @set "class=Bowtie2IndexSet"
# @RdocMethod buildTopHat2TranscriptomeIndexSet
# @alias buildTopHat2TranscriptomeIndexSet
#
# @title "Calls TopHat to build a transcriptome index"
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{gtf}{GtfDataFile to be indexed}
#   \item{outPath}{(optional) Output directory for index and log file.}
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
                                                                             ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'gtf':
  pathnameGTF <- getPathname(gtf)
  pathnameGTF <- Arguments$getReadablePathname(pathnameGTF)

  # Argument 'outPath':
  if (is.null(outPath)) outPath <- file.path(getPath(gtf), "tophat2")
  outPath <- Arguments$getWritablePath(outPath)

  # Argument 'skip':
  skip <- Arguments$getLogical(skip)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose), add=TRUE)
  }


  verbose && enter(verbose, "Building TopHat transcriptome index")

  # The (path) prefix used for the transcripome index set
  prefixName <- getFullName(gtf)
  prefix <- file.path(outPath, prefixName)

  # Check for an existing transcriptome index set
  is <- tryCatch({
    Bowtie2IndexSet$byPrefix(prefix)
  }, error=function(ex) Bowtie2IndexSet())

  # Nothing todo?
  if (skip && isComplete(is)) {
    verbose && cat(verbose, "Transcriptome indexing already done. Skipping.")
    verbose && exit(verbose)
    return(is)
  }

  # Assert TopHat and that it is of a sufficient version
  stopifnot(isCapableOf(aroma.seq, "tophat2"))
  verT <- attr(findTopHat2(), "version")
  if (verT < '2.0.10') throw("TopHat version >= 2.0.10 required")

  # Call TopHat executable
  # (Pre-existing) index for the reference genome
  res <- tophat(getIndexPrefix(this), gtf=gtfFile, outPath=outPath,
                optionsVec=c("--transcriptome-index"=prefixName), ...)

  # Locate index set to return
  is <- Bowtie2IndexSet$byPrefix(prefix)
  verbose && print(verbose, is)

  verbose && exit(verbose)

  is
}) # buildTopHat2TranscriptomeIndexSet()



############################################################################
# HISTORY:
# 2014-07-22 [HB]
# o CLEANUP: Tidied up code and harmonized with the rest of the package.
# 2014-02-04 [TT]
# o Created.
############################################################################

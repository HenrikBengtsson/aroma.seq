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
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasCommas <- function(pathnames, ...) {
    (regexpr(",", pathnames, fixed=TRUE) != -1L);
  } # hasCommas()

  assertNoCommas <- function(pathnames, ...) {
    stopifnot(!any(hasCommas(pathnames)));
  } # assertNoCommas()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'gtf':
  pathnameGTF <- getPathname(gtf)
  pathnameGTF <- Arguments$getReadablePathname(pathnameGTF)
  assertNoCommas(getFullName(gtf))

  # Argument 'outPath':
  if (is.null(outPath)) outPath <- file.path(getParent(getPath(is)), "tophat2")
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
  verbose && print(verbose, "Bowtie2 index set:")
  verbose && print(verbose, this)

  verbose && cat(verbose, "GTF:")
  verbose && print(verbose, gtf)

  # The (path) prefix used for the transcripome index set
  prefixName <- getFullName(gtf)
  prefix <- file.path(outPath, prefixName)
  verbose && cat(verbose, "TopHat2 transcriptome index set prefix: ", prefix)

  # Check for an existing transcriptome and Bowtie2 index set
  # cf. http://ccb.jhu.edu/software/tophat/manual.shtml#tbuild
  tis <- tryCatch({
    Bowtie2IndexSet$byPrefix(prefix)
  }, error=function(ex) Bowtie2IndexSet())

  # Nothing todo?
  if (skip && isComplete(tis)) {
    verbose && cat(verbose, "Transcriptome indexing already done. Skipping.")
    verbose && print(verbose, tis)

    # Assert compatibility
    isCompatibleWith(this, tis, mustWork=TRUE, verbose=less(verbose, 50))

    verbose && exit(verbose)

    return(tis)
  }

  # Assert compatibility
  isCompatibleWith(this, gtf, mustWork=TRUE, verbose=less(verbose, 50))

  # Assert TopHat and that it is of a sufficient version
  stopifnot(isCapableOf(aroma.seq, "tophat2"))
  verT <- attr(findTopHat2(), "version")
  if (verT < '2.0.10') throw("TopHat version >= 2.0.10 required")

  # Call TopHat executable
  # (Pre-existing) index for the reference genome
  res <- tophat(getIndexPrefix(this), gtf=pathnameGTF, outPath=outPath,
                optionsVec=c("--transcriptome-index"="."),
                ..., verbose=less(verbose, 10))
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status: ", status);

  # Locate index set to return
  tis <- Bowtie2IndexSet$byPrefix(prefix)
  verbose && print(verbose, tis)

  # Assert compatibility
  isCompatibleWith(this, tis, mustWork=TRUE, verbose=less(verbose, 50))

  verbose && exit(verbose)

  tis
}) # buildTopHat2TranscriptomeIndexSet()



############################################################################
# HISTORY:
# o ROBUSTNESS: Now buildTopHat2TranscriptomeIndexSet() asserts that the
#   returned index set is compatible with the input index set file.
# 2014-07-22 [HB]
# o ROBUSTNESS: Now buildTopHat2TranscriptomeIndexSet() asserts that the
#   Bowtie2 index set and GTF have compatible sequence names and that
#   the GTF filename has no commas.
# o CLEANUP: Tidied up code and harmonized with the rest of the package.
# 2014-02-04 [TT]
# o Created.
############################################################################

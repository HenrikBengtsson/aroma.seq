###########################################################################/**
# @set "class=Bowtie2IndexSet"
# @RdocMethod buildTopHat2TranscriptomeIndexSet
# @alias buildTopHat2TranscriptomeIndexSet
#
# @title "Calls the TopHat executable to build a transcriptome index on a gene model in a reference genome"
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{tOutDir}{(optional) Output dir for tophat run itself (log files)}
#   \item{gtf}{GtfDataFile to be indexed}
#   \item{tiOutPrefix}{(optional) Output path and prefix for transcriptome index}
#   \item{...}{(Not used)}
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
                                                                             tOutDir=NULL,
                                                                             tiOutPrefix=NULL,
                                                                             ..., verbose=FALSE) {


  organism <- getOrganism(this)

  # Argument gtf
  gtfFile <- getPathname(gtf)
  gtfFile <- Arguments$getReadablePathname(gtfFile)

  # Argument tiOutPrefix
  if (is.null(tiOutPrefix)) {
    tiOutDir <- file.path("annotationData", "organisms", organism, "tophat2")
    tiOutPrefix <- file.path(tiOutDir, sub("[.]g[t|f]f$", "", basename(gtfFile)))
  } else {
    tiOutDir <- dirname(tiOutPrefix)
  }
  Arguments$getWritablePath(tiOutDir)
  stopifnot(!is.null(basename(tiOutPrefix)))

  # Argument tOutDir
  if (is.null(tOutDir)) {
    tOutDir <- file.path(tiOutDir, "tophat2Data")
  }
  Arguments$getWritablePath(tOutDir)

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

  refIdxPrefix <- getIndexPrefix(this)

  # Call TopHat executable

  verbose && enter(verbose, "Building Bowtie2 index set");

  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  verbose && enter(verbose, "Calling tophat()");
  # ( Cmd line call:
  #  'tophat2 -o tOutDir -G gtf --transcriptome-index=tiOutPrefix refIdx' )
  res <- tophat(refIdxPrefix, gtf=gtfFile, optionsVec=c("-o"=tOutDir, "--transcriptome-index"=tiOutPrefix))
  # (0125 TAT - This may not work, because -o is specified in tophat() itself??)

  # Locate existing index files
  res <- tryCatch({
    Bowtie2IndexSet$byPrefix(tiOutPrefix);
  }, error=function(ex) Bowtie2IndexSet());

  verbose && exit(verbose);

  res
}) # buildTopHat2TranscriptomeIndexSet()

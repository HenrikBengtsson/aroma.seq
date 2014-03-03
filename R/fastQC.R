###########################################################################/**
# @RdocDefault fastQC
#
# @title "Calls the FastQC executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnames}{Zero or more input pathnames.}
#   \item{...}{Additional arguments passed to @see "systemFastQC".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# \references{
#  [1] Simon Andrews,
#      FastQC - A quality control tool for high throughput sequence data,
#      March 2014.
#      \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
#  [2] Simon Andrews,
#      FastQC: A quality control application for FastQ data [ANNOUNCEMENT],
#      2010-04-26.
#      \url{http://seqanswers.com/forums/showthread.php?t=4846}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fastQC", "default", function(..., pathnames=character(0L), outPath=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnames':
  pathnames <- lapply(pathnames, FUN=Arguments$getReadablePathname);
  pathnames <- unlist(pathnames, use.names=TRUE);

  # Argument 'outPath':
  if (!is.null(outPath)) {
    outPath <- Arguments$getWritablePath(outPath);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running FastQC");
  verbose && cat(verbose, "Input pathnames:");
  verbose && print(verbose, pathnames);

  args <- list(...);
  if (length(pathnames) > 0L) {
    args <- c(list(pathnames), args);
  }
  if (!is.null(outPath)) {
      str(outPath);
    args <- c(args, sprintf('--outdir "%s"', outPath));
  }
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);
  args <- c(args,verbose=less(verbose, 10));
  res <- do.call(systemFastQC, args=args);

  verbose && exit(verbose);

  res;
}) # fastQC()


############################################################################
# HISTORY:
# 2014-02-28
# o Created from bwaSamse().
############################################################################

###########################################################################/**
# @RdocDefault samtoolsSort
#
# @title "Calls the samtools 'sort' command"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{A SAM/BAM file.}
#   \item{pathnameD}{The destination pathname.}
#   \item{...}{Additional arguments specifying samtools 'sort' switches
#     passed to @see "systemSamtools".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("samtoolsSort", "default", function(pathname, pathnameD, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname);

  # Argument 'pathnameD':
  pathnameD <- Arguments$getWritablePathname(pathnameD);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running samtools 'sort'");

  # Assert that input files are not overwritten
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));

  res <- systemSamtools("sort", ..., shQuote(pathname), shQuote(pathnameD), verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # samtoolsSort()


############################################################################
# HISTORY:
# 2014-03-10 [HB]
# o ROBUSTNESS: Now samtoolsSort() uses shQuote() for all pathnames.
# 2013-11-25
# o Copy from samtoolsView.R
############################################################################


#
# $ samtools sort
#
# Usage:   samtools sort [options] <in.bam> <out.prefix>
#
# Options: -n        sort by read name
#          -f        use <out.prefix> as full file name instead of prefix
#          -o        final output to stdout
#          -l INT    compression level, from 0 to 9 [-1]
#          -@ INT    number of sorting and compression threads [1]
#          -m INT    max memory per thread; suffix K/M/G recognized [768M]
#


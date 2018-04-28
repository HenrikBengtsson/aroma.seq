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
  .stop_if_not(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));

  res <- systemSamtools("sort", ..., shQuote(pathname), shQuote(pathnameD), verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # samtoolsSort()

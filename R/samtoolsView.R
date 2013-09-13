###########################################################################/**
# @RdocDefault samtoolsView
#
# @title "Calls the samtools 'view' command"
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
#   \item{...}{Additional arguments specifying samtools 'view' switches
#     passed to @see "systemSamtools".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("samtoolsView", "default", function(pathname, pathnameD, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Running samtools 'view'");

  # Assert that input files are not overwritten
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));

  res <- systemSamtools("view", "o"=pathnameD, pathname, ..., verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # samtoolsView()


############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################

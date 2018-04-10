###########################################################################/**
# @RdocGeneric convertToBam
# @alias convertToBam.SamDataFile
# @alias convertToBam.SamDataSet
#
# @title "Converts a SAM file (set of files) set into a (sorted and indexed) BAM file (set of files)"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The destination path.}
#  \item{skip}{If @TRUE, already converted files are skipped, otherwise not.}
#  \item{overwrite}{If @TRUE, already converted files are ignored and overwritten.}
#  \item{verbose}{See @see "R.utils::Verbose".}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "BamDataFile" (@see "BamDataSet").
# }
#
# \seealso{
#   Internally @see "Rsamtools::asBam" is utilized.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("convertToBam", "SamDataFile", function(this, path=getPath(this), skip=!overwrite, overwrite=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Converting SAM file to a BAM file");

  pathname <- getPathname(this);
  verbose && cat(verbose, "SAM pathname: ", pathname);

  fullname <- getFullName(this);
  filenameBAM <- sprintf("%s.bam", fullname);
  pathnameBAM <- file.path(path, filenameBAM);
  verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

  # Nothing to do?
  if (skip && isFile(pathnameBAM)) {
    verbose && cat(verbose, "Already converted. Skipping.");
    res <- BamDataFile(pathnameBAM);
    verbose && exit(verbose);
    return(res);
  }

  # Asserts
  stop_if_not(getAbsolutePath(pathnameBAM) != getAbsolutePath(pathname));
  pathnameBAM <- Arguments$getWritablePathname(pathnameBAM, mustNotExist=!overwrite);

  # Converting SAM to BAM, sort and create an index (*.bai)
  verbose && enter(verbose, "Converting using Rsamtools");
  use("Rsamtools")
  pathnameBAMx <- gsub("[.]bam$", "", pathnameBAM);
  verbose && cat(verbose, "BAM destination: ", pathnameBAMx);
  # NB: Rsamtools::asBam() already writes atomically.
  pathnameD <- asBam(pathname, destination=pathnameBAMx,
                     indexDestination=TRUE, overwrite=overwrite);
  verbose && exit(verbose);

  res <- BamDataFile(pathnameBAM);

  verbose && exit(verbose);

  res;
}) # convertToBam() for SamDataFile



setMethodS3("convertToBam", "SamDataSet", function(ds, path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Converting SAM set to a BAM set");

  verbose && cat(verbose, "SAM data set:");
  verbose && print(verbose, ds);
  verbose && cat(verbose, "BAM path: ", path);

  # TO DO: Parallelize. /HB 2013-11-08
  bfList <- list();
  for (ii in seq_along(ds)) {
    df <- ds[[ii]];
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, getName(df), length(ds)));
    pathII <- if (is.null(path)) getPath(df) else path;
    bf <- convertToBam(df, path=pathII, ..., verbose=less(verbose,1));
    bfList[[ii]] <- bf;
    verbose && print(verbose, bf);
    verbose && exit(verbose);
  } # for (ii ...)

  bs <- BamDataSet(bfList);
  bs <- extract(bs, getFullNames(ds), onMissing="error");
  verbose && print(verbose, bs);

  ## TODO: Assert completeness

  verbose && exit(verbose);

  bs;
}) # convertToBam() for SamDataSet

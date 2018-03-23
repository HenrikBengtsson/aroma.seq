###########################################################################/**
# @RdocGeneric sortByPosition
# @alias sortByPosition.BamDataFile
# @alias sortByPosition.BamDataSet
#
# @title "Produce sorted and indexed BAM file(s) from BAM file(s)"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{pathD}{The destination path.}
#  \item{bIndex}{If @TRUE, index file created after sort.}
#  \item{suffix}{Filename suffix for sorted output.}
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
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("sortByPosition", "BamDataFile", function(this, pathD=getPath(this), bIndex=TRUE, suffix=".sorted",
                                                      skip=!overwrite, overwrite=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathD':
  pathD <- Arguments$getWritablePath(pathD);

  # Argument 'bIndex':
  bIndex <- Arguments$getLogical(bIndex)

  # Argument 'suffix'
  stopifnot(is.character(suffix))  ## Should test if single string, and return informative error msg

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

  verbose && enter(verbose, "Sorting BAM file");

  pathname <- getPathname(this);
  verbose && cat(verbose, "Source pathname: ", pathname);

  fullname <- getFullName(this);
  filenameD <- sprintf(paste("%s", suffix, ".bam", sep=""), fullname)
  pathnameD <- file.path(pathD, filenameD);
  verbose && cat(verbose, "Destination pathname: ", pathnameD);

  # Asserts
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=!overwrite);

  # If bam file already sorted, only need to copy .bam and .bai to new destination
  if (skip && isSorted(this)) {
    verbose && cat(verbose, "Already sorted.");
    # Create index if needed
    pathnameI <- paste(pathname, ".bai", sep="")
    if (bIndex && !file.exists(pathnameI)) {
      pathnameI <- indexBam(pathname)
    }
    if (getAbsolutePath(pathnameD) != getAbsolutePath(pathname)) {
      verbose && cat(verbose, "Copying to destination.");
      file.copy(pathname, pathnameD)
      if (bIndex) {
        file.copy(pathnameI, paste(pathnameD, ".bai", sep=""))
      }
    }
    res <- BamDataFile(pathnameD);
    verbose && exit(verbose);
    return(res);
  }

  pathnameDx <- gsub("[.]bam$", "", pathnameD);
  verbose && cat(verbose, "BAM destination: ", pathnameDx);
  # Call Rsamtools::sortBam()
  # - This could hypothetically change pathnameD, but do not check for this (yet)
  # - TODO:  Confirm/modify header to set SORT flag
  pathnameD <- sortBam(file=pathname,
                       destination=pathnameDx)
  if (bIndex) {
    pathnameI <- indexBam(pathnameD)
  }

  res <- BamDataFile(pathnameD);
  verbose && print(verbose, res);
  verbose && exit(verbose);
  res;
}) # sortByPosition() for BamDataFile


setMethodS3("sortByPosition", "BamDataSet", function(ds, path=getPath(ds), suffix=".sorted", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'suffix'
  stopifnot(is.character(suffix))  ## Should test if single string, and return informative error msg

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Sorting BAM data set");

  verbose && cat(verbose, "BAM data set:");
  verbose && print(verbose, ds);
  verbose && cat(verbose, "BAM path: ", path);

  # TO DO: Parallelize. /HB 2013-11-08
  for (ii in seq_along(ds)) {
    df <- ds[[ii]];
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, getName(df), length(ds)));
    bf <- sortByPosition(df, path=path, suffix=suffix, ..., verbose=less(verbose,1));
    verbose && exit(verbose);
  } # for (ii ...)

  bs <- BamDataSet$byPath(path);
  bs <- extract(bs, paste(getFullNames(ds), suffix, sep=""), onMissing="error");
  verbose && print(verbose, bs);

  ## TODO: Assert completeness

  verbose && exit(verbose);

  bs;
}) # sortByPosition() for BamDataSet

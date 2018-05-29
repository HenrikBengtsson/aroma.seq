###########################################################################/**
# @RdocGeneric convertToSam
# @alias convertToSam.BamDataFile
# @alias convertToSam.BamDataSet
#
# @title "Converts a SAM (BAM) file (set of files) set into a BAM (SAM) file (set of files)"
#
# \description{
#   @get "title".
#   When converting to BAM, the generated BAM files are all sorted and indexed as well.
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
#   \code{convertToBam()} returns a (sorted and indexed) @see "BamDataFile" (@see "BamDataSet").
#   \code{convertToSam()} returns a @see "SamDataFile" (@see "SamDataSet").
# }
#
# \seealso{
#   Internally @see "Rsamtools::asBam" or @see "Rsamtools::asSam" is utilized.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("convertToSam", "BamDataFile", function(this, path=getPath(this), skip=!overwrite, overwrite=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path)

  # Argument 'skip':
  skip <- Arguments$getLogical(skip)

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Converting BAM file to a SAM file")

  pathname <- getPathname(this)
  verbose && cat(verbose, "BAM pathname: ", pathname)

  fullname <- getFullName(this)
  filenameSAM <- sprintf("%s.sam", fullname)
  pathnameSAM <- file.path(path, filenameSAM)
  verbose && cat(verbose, "SAM pathname: ", pathnameSAM)

  # Nothing to do?
  if (skip && isFile(pathnameSAM)) {
    verbose && cat(verbose, "Already converted. Skipping.")
    res <- SamDataFile(pathnameSAM)
    verbose && exit(verbose)
    return(res)
  }

  # Asserts
  .stop_if_not(getAbsolutePath(pathnameSAM) != getAbsolutePath(pathname))
  pathnameSAM <- Arguments$getWritablePathname(pathnameSAM, mustNotExist=!overwrite)

  # Converting BAM to SAM
  verbose && enter(verbose, "Converting using Rsamtools")
  asSam <- NULL; rm(list="asSam") # To please R CMD check
  use("Rsamtools")
  pathnameSAMx <- gsub("[.]sam$", "", pathnameSAM)
  verbose && cat(verbose, "SAM destination: ", pathnameSAMx)
  # NB: Rsamtools::asSam() already writes atomically.
  pathnameD <- asSam(pathname, destination=pathnameSAMx,
                     overwrite=overwrite)
  verbose && exit(verbose)

  res <- SamDataFile(pathnameSAM)

  verbose && exit(verbose)

  res
}) # convertToSam() for BamDataFile



setMethodS3("convertToSam", "BamDataSet", function(ds, path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Converting BAM set to a SAM set")

  verbose && cat(verbose, "BAM data set:")
  verbose && print(verbose, ds)
  verbose && cat(verbose, "SAM path: ", path)

  # TO DO: Parallelize. /HB 2013-11-08
  bfList <- list()
  for (ii in seq_along(ds)) {
    df <- ds[[ii]]
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, getName(df), length(ds)))
    pathII <- if (is.null(path)) getPath(df) else path
    bf <- convertToSam(df, path=pathII, ..., verbose=less(verbose,1))
    bfList[[ii]] <- bf
    verbose && print(verbose, bf)
    verbose && exit(verbose)
  } # for (ii ...)

  bs <- SamDataSet(bfList)
  bs <- extract(bs, getFullNames(ds), onMissing="error")
  verbose && print(verbose, bs)

  ## TODO: Assert completeness

  verbose && exit(verbose)

  bs
}) # convertToSam() for BamDataSet

###########################################################################/**
# @RdocGeneric getOutputDataSet
# @alias getOutputDataSet.AromaSeqTransform
# @alias getOutputDataSet.AbstractAlignment
# @alias getOutputDataSet.FastqDownsampler
# @alias getOutputDataSet.PicardDuplicateRemoval
# @alias getOutputDataSet.TopHat2Alignment
# @alias getOutputDataSet.TotalCnBinnedCounting
#
# @title "Gets the (complete or incomplete) processed output data set"
#
# \description{
#   @get "title".
# }
#
# \usage{
#  @usage getOutputDataSet,AromaSeqTransform
#  @usage getOutputDataSet,AbstractAlignment
#  @usage getOutputDataSet,FastqDownsampler
#  @usage getOutputDataSet,PicardDuplicateRemoval
#  @usage getOutputDataSet,TopHat2Alignment
#  @usage getOutputDataSet,TotalCnBinnedCounting
# }
#
# \arguments{
#  \item{onMissing}{A @character string specifying how non-processed files
#   should be returned.
#   If \code{"drop"}, they are ignored and not part of the returned
#   data set.
#   If \code{"dropall"}, @NULL is returned unless all files are processed.
#   If \code{"NA"}, they are represented as a "missing" file.
#   If \code{"error"}, they are not accepted and an exception is thrown.
#  }
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns the output data set containing the same number of files as
#   the input data set, except in the case where argument \code{onMissing}
#   is \code{"drop"} or \code{"dropall"} and one or more files is not
#   processed.
# }
#
# \seealso{
#   This method is utilized by @see "findFilesTodo".
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("getOutputDataSet", "AbstractAlignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path, ...);

  ## Order according to input data set
  ds <- getInputDataSet(this);
  fullnames <- getFullNames(ds);
  bams <- extract(bams, fullnames, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for AbstractAlignment


setMethodS3("getOutputDataSet", "FastqDownsampler", function(this, ...) {
  ## Find all existing output data files
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path, ...);

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  res <- extract(res, fullnames, onMissing="NA");
  res;
}, protected=TRUE) # getOutputDataSet() for FastqDownsampler



setMethodS3("getOutputDataSet", "PicardDuplicateRemoval", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  ds <- getInputDataSet(this);
  path <- getPath(this);
  bams <- byPath(ds, path, ...);

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  bams <- extract(bams, fullnames, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for PicardDuplicateRemoval


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path=path, pattern="accepted_hits.bam$", recursive=TRUE);

  ## Order according to input data set
  # Get the sample names in the found output set
  sampleNames <- getSampleNames(this);
  sampleNamesB <- getPathnames(bams);
  # Workaround: getPathnames() returns NULL, not character(0L), when empty
  if (length(sampleNamesB) == 0L) sampleNamesB <- character(0L);
  sampleNamesB <- basename(dirname(sampleNamesB));
  idxs <- match(sampleNames, sampleNamesB);
  bams <- extract(bams, idxs, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for TopHat2Alignment


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  # AD HOC for now because we are dealing with subdirectories
  # being sample names. /HB 2014-01-10
  paths <- getExpectedOutputPaths(this);
  pathnames <- file.path(paths, "accepted_hits.bam");
  bfList <- lapply(pathnames, FUN=BamDataFile, mustExist=FALSE);
  bams <- BamDataSet(bfList);
  bams <- setFullNamesTranslator(bams, function(names, file, ...) basename(getPath(file)));
  groups  <- getGroups(this);
  fullnames <- names(groups);
  bams <- extract(bams, fullnames, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for TopHat2Alignment


setMethodS3("getOutputDataSet", "HTSeqCounting", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  path <- getPath(this);
  counts <- TabularTextFileSet$byPath(path, ...);
  # FIX ME
  lapply(counts, FUN=function(df) df$readColumnNames <- FALSE);

  ## Order according to input data set
  ds <- getInputDataSet(this);
  fullnames <- getFullNames(ds);
  counts <- extract(counts, fullnames, onMissing=onMissing);

  counts;
}) # getOutputDataSet() for HTSeqCounting


setMethodS3("getOutputDataSet", "TotalCnBinnedCounting", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  # For now, utilize what's already in 'aroma.cn'
  res <- NextMethod("getOutputDataSet", onMissing="drop", .onUnknownArgs="ignore");

  # Don't return NULL
  if (is.null(res)) {
    clazz <- getOutputFileSetClass(this);
    res <- newInstance(clazz, list());
  }

  # Order according to input data set
  ds <- getInputDataSet(this);
  fullnames <- getFullNames(ds);
  res <- extract(res, fullnames, onMissing=onMissing);

  res;
}) # getOutputDataSet() for TotalCnBinnedCounting


############################################################################
# HISTORY:
# 2013-11-15
# o CLEANUP: The different getOutputDataSet() methods no longer have to
#   workaround the special case where the output data set is empty. They
#   also don't have to handle argument 'onMissing'.  All this is now taken
#   care of by extract() for GenericDataFileSet in R.filesets (>= 2.3.3).
# o Extracted all findFilesTodo() methods and document them under the
#   same generic function.
# o Added argument 'onMissing' to getOutputDataSet() for
#   PicardDuplicateRemoval.
# 2013-11-11
# o Added argument 'onMissing' to getOutputDataSet() for AbstractAlignment
#   and TopHat2Alignment.
# 2012-11-26
# o Added argument 'onMissing' to getOutputDataSet() for AbstractAlignment.
# o BUG FIX: getOutputDataSet() for AbstractAlignment would return a data
#   set with "missing" files, if not complete.  Now it only returns the
#   existing files (by default).
############################################################################

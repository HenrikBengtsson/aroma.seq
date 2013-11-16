###########################################################################/**
# @RdocGeneric getOutputDataSet
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
# @synopsis
#
# \arguments{
#  \item{onMissing}{A @character string specifying how non-processed files
#   should be returned.
#   If \code{"drop"}, they are ignored and not part of the returned
#   data set.
#   If \code{"NA"}, they are represented as a "missing" file.
#   If \code{"error"}, they are not accepted and an exception is thrown.
#  }
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns the output data set containing the same number of files as
#   the input data set, except in the case where \code{onMissing="drop"}.
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


  ds <- getInputDataSet(this);

  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path, ...);

  # Special case
  if (length(bams) == 0L) {
    bam <- BamDataFile(NA_character_, mustExist=FALSE);
    bams <- newInstance(bams, list(bam));
  }

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  bams <- extract(bams, fullnames, onMissing="NA");

  # Sanity check
  stopifnot(length(bams) == length(ds));

  exists <- which(unlist(sapply(bams, FUN=isFile)));
  if (length(exists) < length(ds)) {
    if (onMissing == "error") {
      throw("Number of entries in output data set does not match input data set: ", length(exists), " != ", length(ds));
    } else if (onMissing == "drop") {
      bams <- extract(bams, exists);
    }
  }

  # Sanity check
  stopifnot(length(bams) <= length(ds));

  bams;
}) # getOutputDataSet() for AbstractAlignment


setMethodS3("getOutputDataSet", "FastqDownsampler", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path, ...);
  names <- getFullNames(ds);
  idxs <- indexOf(res, names);
  res <- extract(res, idxs, onMissing="NA");
  res;
}, protected=TRUE) # getOutputDataSet() for FastqDownsampler



setMethodS3("getOutputDataSet", "PicardDuplicateRemoval", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  ds <- getInputDataSet(this);

  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path, ...);

  # Special case
  if (length(bams) == 0L) {
    bam <- BamDataFile(NA_character_, mustExist=FALSE);
    bams <- newInstance(bams, list(bam));
  }

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  bams <- extract(bams, fullnames, onMissing="NA");

  # Sanity check
  stopifnot(length(bams) == length(ds));

  exists <- which(unlist(sapply(bams, FUN=isFile)));
  if (length(exists) < length(ds)) {
    if (onMissing == "error") {
      throw("Number of entries in output data set does not match input data set: ", length(exists), " != ", length(ds));
    } else if (onMissing == "drop") {
      bams <- extract(bams, exists);
    }
  }

  # Sanity check
  stopifnot(length(bams) <= length(ds));

  bams;
}) # getOutputDataSet() for PicardDuplicateRemoval


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ds <- getInputDataSet(this);

  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path=path, pattern="accepted_hits.bam$", recursive=TRUE);

  # Special case
  if (length(bams) == 0L) {
    bam <- BamDataFile(NA_character_, mustExist=FALSE);
    bams <- newInstance(bams, list(bam));
  }

  ## Order according to input data set
  # Get the sample names in the found output set
  sampleNamesExpected <- getSampleNames(this);
  sampleNames <- getPathnames(bams);
  sampleNames <- basename(dirname(sampleNames));
  idxs <- match(sampleNamesExpected, sampleNames);
  bams <- extract(bams, idxs, onMissing="NA");

  # Sanity check
  stopifnot(length(bams) == length(ds));
  stopifnot(length(bams) == length(sampleNamesExpected));

  exists <- which(unlist(sapply(bams, FUN=isFile)));
  if (length(exists) < length(ds)) {
    if (onMissing == "error") {
      throw("Number of entries in output data set does not match input data set: ", length(exists), " != ", length(ds));
    } else if (onMissing == "drop") {
      bams <- extract(bams, exists);
    }
  }

  # Sanity check
  stopifnot(length(bams) <= length(ds));

  bams;
}) # getOutputDataSet() for TopHat2Alignment



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

  # Nothing more todo?
  if (onMissing == "drop") {
    return(res);
  }

  ds <- getInputDataSet(this);
  if (length(ds) == 0L) {
    return(res);
  }

  # Special case
  if (length(res) == 0L) {
    className <- getFileClass(res);
    clazz <- Class$forName(className);
    file <- newInstance(clazz, NA_character_, mustExist=FALSE);
    res <- newInstance(res, list(file));
  }

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  res <- extract(res, fullnames, onMissing="NA");

  # Sanity check
  stopifnot(length(res) == length(ds));

  exists <- which(unlist(sapply(res, FUN=isFile)));
  if (length(exists) < length(ds)) {
    if (onMissing == "error") {
      throw("Number of entries in output data set does not match input data set: ", length(exists), " != ", length(ds));
    } else if (onMissing == "drop") {
      res <- extract(res, exists);
    }
  }

  # Sanity check
  stopifnot(length(res) <= length(ds));

  res;
}) # getOutputDataSet() for TotalCnBinnedCounting


############################################################################
# HISTORY:
# 2013-11-15
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

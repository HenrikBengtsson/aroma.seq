###########################################################################/**
# @RdocGeneric getOutputDataSet
# @alias getOutputDataSet.AromaSeqTransform
# @alias getOutputDataSet.AbstractAlignment
# @alias getOutputDataSet.BamDownsampler
# @alias getOutputDataSet.Bowtie2Alignment
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
#  @usage getOutputDataSet,Bowtie2Alignment
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

  ## Sanity check
  if (length(bams) > length(this)) {
    throw(sprintf("Identified more output (%d) than input files (%d): %s",
          length(bams), length(this), sQuote(getPath(bams))));
  }

  ## Assert compatibility with index set
  is <- getIndexSet(this)
  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    bamExists <- isFile(bam)
    if (onMissing != "NA") stopifnot(bamExists)
    if (bamExists) isCompatibleWith(bam, is)
  }

  bams
}, protected=TRUE)


setMethodS3("getOutputDataSet", "BamDownsampler", function(this, ...) {
  ## Find all existing output data files
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path, ...);

  ## Order according to input data set
  fullnames <- getFullNames(ds);
  res <- extract(res, fullnames, onMissing="NA");
  res;
}, protected=TRUE) # getOutputDataSet() for BamDownsampler


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


setMethodS3("getOutputDataSet", "Bowtie2Alignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);

  path <- getPath(this);
  filenames <- sprintf("%s.bam", getSampleNames(this))
  pathnames <- file.path(path, filenames);
  bfList <- lapply(pathnames, FUN=BamDataFile, mustExist=FALSE);
  bams <- BamDataSet(bfList);

  groups  <- getGroups(this);
  fullnames <- names(groups);
  bams <- extract(bams, fullnames, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for Bowtie2Alignment


setMethodS3("getOutputDataSet", "HTSeqCounting", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  path <- getPath(this);
  counts <- HTSeqCountDataSet$byPath(path, ...);

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

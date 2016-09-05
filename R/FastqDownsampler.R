###########################################################################/**
# @RdocClass FastqDownsampler
#
# @title "The FastqDownsampler class"
#
# \description{
#  @classhierarchy
#
#  ...
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "FastqDataSet".}
#  \item{subset}{An @integer specifying the total number of reads to sample,
#    or a @double specifying the fraction of total number of reads to sample.}
#  \item{seed}{An (optional) @integer specifying the random seed to be
#     set before sampling indices.  The random seed is set to its original
#     state when exiting.  If @NULL, it is not set.}
#  \item{...}{Additional arguments passed to @see "AromaSeqTransform".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#  Internally, the @see "ShortRead::FastqSampler" method is used.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FastqDownsampler", function(dataSet=NULL, subset=1e6, seed=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "FastqDataSet");

    # Argument 'subset':
    if (length(subset) == 1L) {
      subset <- Arguments$getNumeric(subset, range=c(0,Inf));
      if (subset <= 1) {
        subset <- Arguments$getDouble(subset, range=c(0,1));
      } else {
        subset <- Arguments$getInteger(subset, range=c(1,Inf));
      }
    } else {
      throw("Not yet implemented.");
      subset <- Arguments$getIndex(subset);
    }

    # Argument 'seed':
    if (!is.null(seed)) seed <- Arguments$getInteger(seed);
  } # if (!is.null(dataSet))

  extend(AromaSeqTransform(dataSet=dataSet, subset=subset, seed=seed, ...), "FastqDownsampler");
})


setMethodS3("getSampleSize", "FastqDownsampler", function(this, df, ...) {
  params <- getParameters(this);
  subset <- params$subset;
  if (subset <= 1) {
    n <- subset * nbrOfSeqs(df);
    n <- Arguments$getInteger(n);
  } else {
    n <- subset;
  }
  n;
}, protected=TRUE);


setMethodS3("getAsteriskTags", "FastqDownsampler", function(this, ...) {
  params <- getParameters(this);
  sprintf("n=%g", params$subset);
}, protected=TRUE);


setMethodS3("getRootPath", "FastqDownsampler", function(this, ...) {
  "fastqData";
}, protected=TRUE);


setMethodS3("process", "FastqDownsampler", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Downsampling FASTQ data set");
  ds <- getInputDataSet(this);
  verbose && print(verbose, ds);

  params <- getParameters(this);
  seed <- params$seed;
  verbose && cat(verbose, "Random seed: ", seed);

  if (isPaired(ds)) {
    throw(sprintf("%s does not yet support paired-end FASTQ data sets: %s",
                  class(this)[1L], getPath(ds)));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this)
  skip <- !force

  res <- listenv()
  for (ii in seq_along(ds)) {
    df <- ds[[ii]]
    fullname <- getFullName(df)

    verbose && enter(verbose, sprintf("FASTQ downsampling sample #%d ('%s') of %d", ii, fullname, length(ds)))

    ext <- tools::file_ext(getFilename(df))
    filename <- sprintf("%s.%s", fullname, ext)
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE)
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      res[[ii]] <- pathname
      verbose && exit(verbose)
    }

    res[[ii]] %<-% {
      verbose && print(verbose, df)

      n <- getSampleSize(this, df)
      verbose && printf(verbose, "Sample size: %d\n", n)

      if (isFile(pathname)) file.remove(pathname)
      dfT <- writeSample(df, n=n, ordered=FALSE, seed=seed, pathname=pathname)
      verbose && print(verbose, dfT)

      getPathname(dfT)
    }

    verbose && exit(verbose)
  } ## for (ii ...)

  ## Resove futures
  res <- resolve(res)

  res <- getOutputDataSet(this, verbose=less(verbose, 1))

  verbose && exit(verbose)

  invisible(res)
})


############################################################################
# HISTORY:
# 2014-04-18
# o BUG FIX: FastqDownsampler did not respect fullname translators.
# 2013-11-16
# o CLEANUP: Dropped several methods now taken care of by super class.
# 2013-09-12
# o Now process() for FastqDownsampler utilizes parallel processing.
# o CLEANUP: Improved previous "mockup" code of FastqDownsampler.
# 2013-09-03
# o ROBUSTNESS: Now process() for FastqDownsampler gives an error
#   if the data set is paired-end; will implement later.
# 2013-07-01
# o Created.
############################################################################

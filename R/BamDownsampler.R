###########################################################################/**
# @RdocClass BamDownsampler
#
# @title "The BamDownsampler class"
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
#  \item{dataSet}{An @see "BamDataSet".}
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
#  Internally, the @see "Rsamtools::BamSampler" method is used.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("BamDownsampler", function(dataSet=NULL, subset=1e6, seed=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");

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

  extend(AromaSeqTransform(dataSet=dataSet, subset=subset, seed=seed, ...), "BamDownsampler");
})


setMethodS3("getSampleSize", "BamDownsampler", function(this, df, ...) {
  params <- getParameters(this);
  subset <- params$subset;
  if (subset <= 1) {
    n <- subset * nbrOfReads(df);
    n <- Arguments$getInteger(n);
  } else {
    n <- subset;
  }
  n;
}, protected=TRUE);


setMethodS3("getAsteriskTags", "BamDownsampler", function(this, ...) {
  params <- getParameters(this);
  sprintf("n=%g", params$subset);
}, protected=TRUE);


setMethodS3("getRootPath", "BamDownsampler", function(this, ...) {
  "bamData";
}, protected=TRUE);


setMethodS3("process", "BamDownsampler", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Downsampling BAM data set")
  ds <- getInputDataSet(this)
  verbose && print(verbose, ds)

  params <- getParameters(this)
  seed <- params$seed
  verbose && cat(verbose, "Random seed: ", seed)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the BAM files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  skip <- !force

  res <- listenv()
  for (ii in seq_along(ds)) {
    df <- ds[[ii]]
    fullname <- getFullName(df)

    verbose && enter(verbose, sprintf("Downsampling BAM #%d ('%s') of %d", ii, fullname, length(ds)))

    ext <- tools::file_ext(getFilename(df))
    filename <- sprintf("%s.%s", fullname, ext)
    pathname <- Arguments$getWritablePathname(filename, path=getPath(this), mustNotExist=FALSE)

    ## Already processed?
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      verbose && exit(verbose)
      res[[ii]] <- NA
      next
    }

    res[[ii]] %<-% {
      verbose && print(verbose, df)
      n <- getSampleSize(this, df)
      verbose && printf(verbose, "Sample size: %d\n", n)
      if (isFile(pathname)) file.remove(pathname)
      dfT <- writeSample(df, n=n, seed=seed, pathname=pathname, verbose=verbose)
      verbose && print(verbose, dfT)
      dfT
    } ## %<-%

    verbose && exit(verbose)
  } ## for (ii ...)

  ds <- NULL  ## Not needed anymore

  ## Resolve all futures
  res <- resolve(res)

  res <- getOutputDataSet(this, verbose=less(verbose, 1))

  verbose && exit(verbose)

  invisible(res)
})

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
setConstructorS3("BamDownsampler", function(dataSet=NULL, subset=1e6, ...) {
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
  } # if (!is.null(dataSet))

  extend(AromaSeqTransform(dataSet=dataSet, subset=subset, ...), "BamDownsampler");
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
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Downsampling BAM data set");

  ds <- getInputDataSet(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the BAM files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(df, sampler, path, ..., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

    # Argument 'sampler':
    stopifnot(is.function(sampler));

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }


    verbose && enter(verbose, "BAM downsampling of one sample");

    fullname <- getFullName(df);
    ext <- getFilenameExtension(df);
    filename <- sprintf("%s.%s", fullname, ext);
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.");
      verbose && exit(verbose);
      return(invisible(list(pathnameFQ=pathname, n=NA_integer_)));
    }

    verbose && print(verbose, df);
    n <- sampler(df);
    verbose && printf(verbose, "Sample size: %d\n", n);
    if (isFile(pathname)) {
      file.remove(pathname);
    }

    dfT <- writeSample(df, n=n, pathname=pathname);
    verbose && print(verbose, dfT);

    # Not needed anymore
    df <- dfT <- NULL;

    verbose && exit(verbose);

    invisible(list(pathnameFQ=pathname, n=n));
  }, sampler=function(df) { getSampleSize(this, df) }, path=getPath(this), skip=!force, verbose=verbose) # dsApply()

  # Not needed anymore
  ds <- NULL;

  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2014-04-18
# o Created from FastqDownsampler.R.
############################################################################

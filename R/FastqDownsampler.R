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
#  \item{tags}{Tags for the output data set.}
#  \item{...}{Not used.}
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
setConstructorS3("FastqDownsampler", function(dataSet=NULL, subset=1e6, tags="*", ...) {
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
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(...);

  # TODO: Turn into one of the Transformer classes
  extend(Object(), c("FastqDownsampler"),
    .ds = dataSet,
    .subset = subset,
    .tags = tags
  );
})


setMethodS3("as.character", "FastqDownsampler", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  ds <- getInputDataSet(this);
  s <- c(s, "Input data set:");
  s <- c(s, as.character(ds));

  s <- c(s, sprintf("Output path: %s", getPath(this)));
  s <- c(s, sprintf("Added tags: %s", getTags(this, collapse=TRUE)));

  # Algorithm parameters
  s <- c(s, sprintf("Sample size: %s", getSampleSize(this)));

  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)



setMethodS3("getSampleSize", "FastqDownsampler", function(this, df, ...) {
  subset <- this$.subset;
  if (subset <= 1) {
    n <- subset * nbrOfSeqs(df);
    n <- Arguments$getInteger(n);
  } else {
    n <- subset;
  }
  n;
}, protected=TRUE);


setMethodS3("getInputDataSet", "FastqDownsampler", function(this, ...) {
  this$.ds;
}, protected=TRUE)


setMethodS3("getOutputRootPath", "FastqDownsampler", function(this, ...) {
  "fastqData";
}, protected=TRUE);


setMethodS3("getAsteriskTags", "FastqDownsampler", function(this, ...) {
  sprintf("n=%g", this$.subset);
}, protected=TRUE);


setMethodS3("getTags", "FastqDownsampler", function(this, collapse=FALSE, ...) {
  tags <- this$.tags;
  if (is.null(tags)) return(tags);
  tags <- unlist(strsplit(tags, split=",", fixed=TRUE));
  tags[tags == "*"] <- getAsteriskTags(this);
  if (collapse) {
    tags <- paste(tags, collapse=",");
  }
  tags;
}, protected=TRUE);


setMethodS3("getPath", "FastqDownsampler", function(this, ...) {
  rootPath <- getOutputRootPath(this);
  ds <- getInputDataSet(this);

  # Output dataset name
  fullname <- getFullName(ds);
  tags <- getTags(this, collapse=TRUE);
  fullname <- paste(c(fullname, tags), collapse=",");

  chipType <- basename(getPath(ds));
  path <- file.path(rootPath, fullname, chipType);

  path <- Arguments$getWritablePath(path);
  path;
}, protected=TRUE)



setMethodS3("getOutputDataSet", "FastqDownsampler", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path, ...);
  names <- getFullNames(ds);
  idxs <- indexOf(res, names);
  res <- extract(res, idxs, onMissing="NA");
  res;
}, protected=TRUE)



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

  if (isPaired(ds)) {
    throw(sprintf("%s does not yet support paired-end FASTQ data sets: %s",
                  class(this)[1L], getPathname(ds)));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(df, sampler, path, ..., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "FastqDataFile");

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


    verbose && enter(verbose, "FASTQ downsampling of one sample");

    filename <- getFilename(df);
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

    dfT <- writeSample(df, n=n, ordered=FALSE, pathname=pathname);
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
# 2013-09-12
# o Now process() for FastqDownsampler utilizes dsApply().
# o CLEANUP: Improved previous "mockup" code of FastqDownsampler.
# 2013-09-03
# o ROBUSTNESS: Now process() for FastqDownsampler gives an error
#   if the data set is paired-end; will implement later.
# 2013-07-01
# o Created.
############################################################################

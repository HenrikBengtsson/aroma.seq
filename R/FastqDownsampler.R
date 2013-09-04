setConstructorS3("FastqDownsampler", function(dataSet=NULL, subset=1e6, ...) {
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
    .subset = subset
  );
})


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

setMethodS3("getOutputPath", "FastqDownsampler", function(this, ...) {
  ds <- getInputDataSet(this);
  fullname <- getFullName(ds);
  nTag <- sprintf("n=%g", this$.subset);
  fullname <- paste(c(fullname, nTag), collapse=",");
  rootPath <- getOutputRootPath(this);
  path <- file.path(rootPath, fullname);
  path <- Arguments$getWritablePath(path);
  path;
}, protected=TRUE)


setMethodS3("process", "FastqDownsampler", function(this, ..., force=FALSE) {
  ds <- getInputDataSet(this);
  path <- getOutputPath(this);

  if (isPaired(ds)) {
    throw(sprintf("%s does not yet support paired-end FASTQ data sets: %s",
                  class(this)[1L], getPathname(ds)));
  }

  dfTList <- list();
  for (ii in seq_along(this)) {
    df <- getFile(ds, ii);
    filename <- getFilename(df);
    pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=FALSE);
    if (!force && isFile(pathname)) {
      dfT <- newInstance(df, pathname);
      dfTList[[ii]] <- dfT;
      next;
    }
    n <- getSampleSize(this, df);
    print(n);
    if (isFile(pathname)) {
      file.remove(pathname);
    }
    dfT <- writeSample(df, n=n, ordered=FALSE, pathname=pathname);
    print(dfT);
    dfTList[[ii]] <- dfT;

    # Not needed anymore
    df <- dfT <- NULL;
  } # for (ii ...)

  res <- newInstance(ds, dfTList);

  # Not needed anymore
  ds <- dfTList <- NULL;

  res;
})


############################################################################
# HISTORY:
# 2013-09-03
# o ROBUSTNESS: Now process() for FastqDownsampler gives an error
#   if the data set is paired-end; will implement later.
# 2013-07-01
# o Created.
############################################################################

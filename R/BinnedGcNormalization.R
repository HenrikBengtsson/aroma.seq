###########################################################################/**
# @RdocClass BinnedGcNormalization
#
# @title "The abstract BinnedGcNormalization class"
#
# \description{
#  @classhierarchy
#
# }
# 
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "aroma.core::AromaUnitTotalCnBinarySet".}
#  \item{...}{Arguments passed to @see "aroma.core::AromaTransform".}
#  \item{.reqSetClass}{(internal only)}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 
setConstructorS3("BinnedGcNormalization", function(dataSet=NULL, ..., .reqSetClass="AromaUnitTotalCnBinarySet") {
  if (!is.null(dataSet)) {
  }

  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "BinnedGcNormalization"
  );
}, abstract=TRUE)


setMethodS3("as.character", "BinnedGcNormalization", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  clazz <- class(s);
  unc <- getGcContentFile(this);
##  s <- c(s, "GcContentFile:");
  s <- c(s, as.character(unc));
  class(s) <- clazz;
  s;
}, protected=TRUE)



setMethodS3("getParameters", "BinnedGcNormalization", function(this, ...) {
  params <- NextMethod("getParameters");
  params;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "BinnedGcNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class-specific tags

  params <- getParameters(this);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } 

  tags;
}, protected=TRUE) 



setMethodS3("getRootPath", "BinnedGcNormalization", function(this, ...) {
  "smoothCnData";
}, protected=TRUE)



setMethodS3("getPath", "BinnedGcNormalization", function(this, create=TRUE, ...) {
  path <- NextMethod("getPath", create=FALSE);
  parent <- dirname(path);

  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(parent, chipType);

  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)


setMethodS3("getGcContentFile", "BinnedGcNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  ugp <- getAromaUgpFile(ds);
  unc <- getAromaUncFile(ugp);
  unc;
})


setMethodS3("getOutputFileExtension", "BinnedGcNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  df <- getFile(ds, 1);
  ext <- getFilenameExtension(df);
  sprintf(".%s", ext);
}, protected=TRUE)


setMethodS3("getOutputFileSetClass", "BinnedGcNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  className <- class(ds)[1];
  Class$forName(className);
}, protected=TRUE)


setMethodS3("getOutputFileClass", "BinnedGcNormalization", function(this, ...) {
  setClass <- getOutputFileSetClass(this, ...);
  setInstance <- newInstance(setClass);
  className <- getFileClass(setInstance);
  Class$forName(className);
}, protected=TRUE)


setMethodS3("getOutputDataSet0", "BinnedGcNormalization", function(this, pattern=NULL, className=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern=pattern);
  }


  verbose && enter(verbose, "Retrieving existing set of output files");
  ds <- getInputDataSet(this);
  outPath <- getPath(this);
  if (is.null(className)) {
    clazz <- getOutputFileSetClass(this);
    className <- getName(clazz);
  }
  verbose && cat(verbose, "Class: ", className);

  path <- getPath(this);
  verbose && cat(verbose, "Path: ", path);

  if (is.null(pattern)) {
    # Default filename pattern find non-private (no dot prefix) files with
    # the same file name extension as the input data set.
    fileExt <- getOutputFileExtension(this);
    fileExt <- c(fileExt, tolower(fileExt), toupper(fileExt));
    fileExt <- sprintf("(%s)", paste(unique(fileExt), collapse="|"));
    verbose && cat(verbose, "Expected file extensions: ", fileExt);
    pattern <- sprintf("^[^.].*%s$", fileExt);
  }
  verbose && cat(verbose, "Pattern: ", pattern);

  verbose && enter(verbose, sprintf("Calling %s$forName()", className));
  clazz <- Class$forName(className);
  args <- list(path=path, pattern=pattern, ...);
  verbose && str(verbose, args);
  args$verbose <- less(verbose);
  staticMethod <- clazz$byPath;
  dsOut <- do.call("staticMethod", args=args);
  rm(staticMethod, args); # Not needed anymore
  verbose && exit(verbose);

  verbose && exit(verbose);

  dsOut;
}, protected=TRUE)



setMethodS3("process", "BinnedGcNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (isDone(this)) {
    dsOut <- getOutputDataSet(this);
    return(invisible(dsOut));
  }

  verbose && enter(verbose, "Normalizing for binned GC content effects");

  params <- getParameters(this);

  verbose && print(verbose, "Input data set:");
  ds <- getInputDataSet(this);
  verbose && print(verbose, ds);

  verbose && enter(verbose, "Extracting GC content per bin");
  unc <- getGcContentFile(this);
  verbose && print(verbose, unc);

  gc <- getGcContent(unc);
  verbose && str(verbose, gc);
  verbose && summary(verbose, gc);
  verbose && exit(verbose);

  nbrOfUnits <- nbrOfUnits(unc);

  # Get Class object for the output files
  clazz <- getOutputFileClass(this);

  # Get the filename extension for output files
  ext <- getOutputFileExtension(this);


  for (ii in seq_along(ds)) {
    df <- getFile(ds, ii);
    name <- getFullName(df);
    verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", ii, name, length(ds)));

    path <- getPath(this);
    fullname <- getFullName(df);
    filename <- sprintf("%s%s", fullname, ext);
    pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                         mustExist=FALSE);
    verbose && cat(verbose, "Output pathname: ", pathname);

    if (isFile(pathname)) {
      dfOut <- newInstance(clazz, filename=pathname);
      if (nbrOfUnits != nbrOfUnits(dfOut)) {
        throw("The number of units in existing output file does not match the number of units in the output file: ", nbrOfUnits, " != ", nbrOfUnits(dfOut));
      }
      verbose && cat(verbose, "Skipping already existing output file.");
      verbose && exit(verbose);
      next;
    }

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    verbose && print(verbose, df);

    y <- extractMatrix(df, drop=TRUE, verbose=less(verbose, 10));
    verbose && cat(verbose, "Signals:");
    verbose && str(verbose, y);

    verbose && enter(verbose, "Normalizing signals (on the log scale) for GC content");
    ly <- log2(y);
    targetFcn <- function(...) 1;
    lyN <- normalizeGcContent(ly, gcContent=gc, targetFcn=targetFcn, .isLogged=TRUE, .returnFit=TRUE);
    yN <- 2^lyN;
    fit <- attr(lyN, "modelFit");
    verbose && cat(verbose, "Model fit:");
    verbose && str(verbose, fit);

    rm(y, ly, lyN);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing normalized signals");
    verbose && cat(verbose, "Pathname: ", pathname);

    paramsT <- params;
    footer <- list(
      sourceDataFile=list(
        fullname=getFullName(df), 
        platform=getPlatform(df), 
        chipType=getChipType(df), 
        checksum=getChecksum(df)
      ), parameters=list(
        annotation=list(
          fullname=getFullName(unc),
          platform=getPlatform(unc),
          chipType=getChipType(unc),
          checksum=getChecksum(unc)
        ),
        params=paramsT
      )
    );

    platform <- getPlatform(df);
    chipType <- getChipType(df);

    dfOut <- clazz$allocate(filename=pathnameT, nbrOfRows=nbrOfUnits, 
                            platform=platform, chipType=chipType, 
                            footer=footer, verbose=less(verbose, 50));

    dfOut[,1] <- yN;
    rm(yN);

    # Renaming temporary file
    pathname <- popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose); # Storing

    verbose && exit(verbose);
  } # for (ii ...)

  verbose && exit(verbose);

  dsOut <- getOutputDataSet(this);
  invisible(dsOut);
})



setMethodS3("getOutputFiles", "BinnedGcNormalization", function(this, ...) {
  NextMethod("getOutputFiles", pattern=".*[.]asb$");
}, protected=TRUE) 



############################################################################
# HISTORY:
# 2012-10-18
# o Created.
############################################################################

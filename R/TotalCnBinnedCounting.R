###########################################################################/**
# @RdocClass TotalCnBinnedCounting
#
# @title "The TotalCnBinnedCounting class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "aroma.cn::TotalCnSmoothing".}
#  \item{.reqSetClass}{(internal) ...}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("TotalCnBinnedCounting", function(..., .reqSetClass="BamDataSet") {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  extend(TotalCnSmoothing(..., .reqSetClass=.reqSetClass), "TotalCnBinnedCounting");
})


setMethodS3("getExpectedOutputFullnames", "TotalCnBinnedCounting", function(this, ...) {
  names <- NextMethod("getExpectedOutputFullnames");
  names <- paste(names, "counts", sep=",");
  names;
}, protected=TRUE)


setMethodS3("getOutputFileSetClass", "TotalCnBinnedCounting", function(this, ...) {
  AromaUnitTotalCnBinarySet;
}, protected=TRUE)

setMethodS3("getOutputFileExtension", "TotalCnBinnedCounting", function(this, ...) {
  ",counts.asb";
}, protected=TRUE)


setMethodS3("smoothRawCopyNumbers", "TotalCnBinnedCounting", function(this, rawCNs, target, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Counting one set of copy numbers");
  verbose && print(verbose, rawCNs);

  x <- getPositions(rawCNs);
  verbose && cat(verbose, "Read positions to be binned:");
  verbose && str(verbose, x);

  verbose && cat(verbose, "Argument 'target':");
  verbose && str(verbose, target);

  # Setting up arguments
  params <- getParameters(this);
  targetUgp <- params$targetUgp;
  params$targetUgp <- NULL;

  xOut <- target$xOut;
  verbose && cat(verbose, "Target loci: ", hpaste(xOut));
  by <- median(diff(sort(xOut)), na.rm=TRUE);

  verbose && cat(verbose, "Distance between target loci: ", by);
  bx <- c(xOut[1]-by/2, xOut+by/2);

  verbose && cat(verbose, "Bins:");
  verbose && str(verbose, bx);

  args <- c(list(), params, list(x=x, bx=bx), ...);

  # Keep only known arguments
  knownArguments <- names(formals(binCounts.default));
  keep <- is.element(names(args), knownArguments);
  args <- args[keep];

  verbose && cat(verbose, "Calling binCounts() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 20);
  yS <- do.call("binCounts", args=args);
  verbose && cat(verbose, "Bin counts:");
  verbose && str(verbose, yS);

  smoothCNs <- RawCopyNumbers(yS, x=xOut);

  verbose && exit(verbose);

  smoothCNs;
}, protected=TRUE)



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
})

setMethodS3("findFilesTodo", "TotalCnBinnedCounting", function(this, ...) {
  res <- getOutputDataSet(this, onMissing="NA");
  isFile <- unlist(sapply(res, FUN=isFile), use.names=FALSE);
  todo <- !isFile;
  todo <- which(todo);
  if (length(todo) > 0L) {
    ds <- getInputDataSet(this);
    names(todo) <- getNames(ds[todo]);
  }
  todo;
})


setMethodS3("process", "TotalCnBinnedCounting", function(this, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Smoothing copy-number towards set of target loci");

  params <- getParameters(this);

  verbose && print(verbose, "Input data set:");
  ds <- getInputDataSet(this);
  verbose && print(verbose, ds);

  verbose && enter(verbose, "Identifying all target positions");
  targetList <- getTargetPositions(this, ...);
  nbrOfChromosomes <- length(targetList);
  verbose && str(verbose, targetList);

  targetUgp <- params$targetUgp;
  platform <- getPlatform(targetUgp);
  chipType <- getChipType(targetUgp);
  nbrOfUnits <- nbrOfUnits(targetUgp);
  # Not needed anymore
  targetUgp <- NULL;
  verbose && cat(verbose, "Total number of target units:", nbrOfUnits);
  verbose && exit(verbose);

  # Get Class object for the output files
  clazz <- getOutputFileClass(this);

  # Get the filename extension for output files
  ext <- getOutputFileExtension(this);


  nbrOfArrays <- length(ds);
  for (kk in seq_along(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array %d ('%s') of %d",
                                            kk, getName(df), nbrOfArrays));

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

    verbose && print(verbose, df);

    # Preallocate vector
    M <- rep(as.double(NA), times=nbrOfUnits);

    verbose && enter(verbose, "Reading and smoothing input data");
    for (cc in seq_along(targetList)) {
      target <- targetList[[cc]];
      chromosome <- target$chromosome;
      chrTag <- sprintf("Chr%02d", chromosome);

      verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d",
                                               cc, chrTag, nbrOfChromosomes));
      verbose && cat(verbose, "Extracting raw CNs:");
      rawCNs <- extractRawCopyNumbers(df, chromosome=chromosome,
                                                  verbose=less(verbose, 10));
      verbose && print(verbose, rawCNs);
      verbose && summary(verbose, rawCNs);

      verbose && cat(verbose, "Smoothing CNs:");
      verbose && cat(verbose, "Target positions:");
      verbose && str(verbose, target$xOut);

      smoothCNs <- smoothRawCopyNumbers(this, rawCNs=rawCNs,
                                        target=target, verbose=verbose);

      verbose && print(verbose, smoothCNs);
      verbose && summary(verbose, smoothCNs);

      M[target$units] <- getSignals(smoothCNs);
      verbose && exit(verbose);
    } # for (cc ...)

    verbose && cat(verbose, "Smoothed CNs across all chromosomes:");
    verbose && str(verbose, M);
    verbose && summary(verbose, M);
    verbose && printf(verbose, "Missing values: %d (%.1f%%) out of %d\n",
                   sum(is.na(M)), 100*sum(is.na(M))/nbrOfUnits, nbrOfUnits);
    verbose && exit(verbose);

    verbose && enter(verbose, "Storing smoothed data");
    verbose && cat(verbose, "Pathname: ", pathname);

    params2 <- params;
    params2[["targetUgp"]] <- NULL;
    footer <- list(
      sourceDataFile=list(
        fullname=getFullName(df),
        platform=getPlatform(df),
        chipType=getChipType(df),
        checksum=getChecksum(df)
      ), parameters=list(
        targetUgp=list(
          fullname=getFullName(params$targetUgp),
          platform=getPlatform(params$targetUgp),
          chipType=getChipType(params$targetUgp),
          checksum=getChecksum(params$targetUgp)
        ),
        params=params2
      )
    );

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    dfOut <- clazz$allocate(filename=pathnameT, nbrOfRows=nbrOfUnits,
                            platform=platform, chipType=chipType,
                            footer=footer, verbose=less(verbose, 50));

    dfOut[,1] <- M;
    # Not needed anymore
    M <- NULL;

    # Renaming temporary file
    pathname <- popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose); # Storing

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  dsOut <- getOutputDataSet(this);
  invisible(dsOut);
}) # process()



setMethodS3("extractRawCopyNumbers", "BamDataFile", function(this, chromosome, ..., verbose=FALSE) {
  require("GenomicRanges") || throw("Package not loaded: GenomicRanges");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting raw \"copy numbers\"");
  verbose && cat(verbose, "Chromosome index: ", chromosome);

  targetLabels <- names(getTargets(this));
  chrLabel <- targetLabels[chromosome];  # AD HOC. /HB 2012-10-11
  verbose && cat(verbose, "Chromosome label: ", chrLabel);

  # Read (start) positions on current chromosome
  gr <- GRanges(seqnames=Rle(chrLabel), ranges=IRanges(-500e6, +500e6));
  x <- readReadPositions(this, which=gr, verbose=less(verbose, 10))$pos;
  verbose && cat(verbose, "Read positions:");
  verbose && str(verbose, x);

  y <- rep(1.0, times=length(x));
  cn <- RawCopyNumbers(y, x=x, chromosome=chromosome);

  verbose && cat(verbose, "Read data:");
  verbose && cat(verbose, cn);

  verbose && exit(verbose);

  cn;
}, protected=TRUE) # extractRawCopyNumbers()


setMethodS3("getPlatform", "BamDataSet", function(this, ...) {
  getPlatform(getOneFile(this));
})

setMethodS3("getChipType", "BamDataSet", function(this, ...) {
  getChipType(getOneFile(this));
})

setMethodS3("getPlatform", "BamDataFile", function(this, ...) {
  "NGS";
})

setMethodS3("getChipType", "BamDataFile", function(this, ...) {
  basename(getPath(this));
})

setMethodS3("getFilenameExtension", "BamDataFile", function(this, ...) {
  "bam";
}, protected=TRUE)


############################################################################
# HISTORY:
# 2012-10-11
# o Created from TotalCnBinnedSmoothing.R.
############################################################################

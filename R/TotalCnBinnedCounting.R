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
  use("aroma.cn")
  TotalCnSmoothing <- aroma.cn::TotalCnSmoothing

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
  ns <- getNamespace("matrixStats");
  binCounts <- get("binCounts", envir=ns, mode="function");
  knownArguments <- names(formals(binCounts));
  rm(list="binCounts");
  keep <- is.element(names(args), knownArguments);
  args <- args[keep];

  verbose && cat(verbose, "Calling binCounts() with arguments:");
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 20);
  yS <- do.call(binCounts, args=args);
  verbose && cat(verbose, "Bin counts:");
  verbose && str(verbose, yS);

  smoothCNs <- RawCopyNumbers(yS, x=xOut);

  verbose && exit(verbose);

  smoothCNs;
}, protected=TRUE)




setMethodS3("process", "TotalCnBinnedCounting", function(this, ..., force=FALSE, verbose=FALSE) {
  use("aroma.cn")
  getTargetPositions <- aroma.cn::getTargetPositions

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Smoothing copy-number towards set of target loci");

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  if (force) {
    todo <- seq_along(ds);
  } else {
    todo <- findFilesTodo(this, verbose=less(verbose, 1));
    # Already done?
    if (length(todo) == 0L) {
      verbose && cat(verbose, "Already done. Skipping.");
      res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
      verbose && exit(verbose);
      return(invisible(res));
    }
  }

  params <- getParameters(this);
  verbose && cat(verbose, "Method parameters:");
  verbose && str(verbose, params);

  verbose && enter(verbose, "Identifying all target positions");
  targetList <- getTargetPositions(this, ...);
  nbrOfChromosomes <- length(targetList);
  verbose && str(verbose, targetList);

  targetUgp <- params$targetUgp;
  platform <- getPlatform(targetUgp);
  chipType <- getChipType(targetUgp);
  nbrOfUnits <- nbrOfUnits(targetUgp);
  params$targetUgp <- NULL;
  parameters <- list(
    targetUgp=list(
      fullname=getFullName(targetUgp),
      platform=getPlatform(targetUgp),
      chipType=getChipType(targetUgp),
      nbrOfUnits=nbrOfUnits(targetUgp),
      checksum=getChecksum(targetUgp)
    ),
    params=params
  )
  # Not needed anymore
  targetUgp <- params <- NULL;
  verbose && cat(verbose, "Total number of target units: ", nbrOfUnits);
  verbose && exit(verbose);

  # Get Class object for the output files
  clazz <- getOutputFileClass(this);

  # Get the filename extension for output files
  ext <- getOutputFileExtension(this);

  # Output directory
  path <- getPath(this)
  verbose && cat(verbose, "Number of items to be processed: ", length(todo))

  ## Number of units
  nbrOfUnits <- parameters$targetUgp$nbrOfUnits

  ## FIXME: For some reason the below does not work with
  ##        parallel futures. /HB 2016-03-12
  oplan <- future::plan("eager")
  on.exit(future::plan(open), add=TRUE)

  res <- listenv()
  for (kk in seq_along(todo)) {
    ii <- todo[kk]
    df <- ds[[ii]]
    fullname <- getFullName(df)

    verbose && enter(verbose, sprintf("Item #%d ('%s') of %d", kk, fullname, length(todo)))

    filename <- sprintf("%s%s", fullname, ext)
    pathname <- Arguments$getReadablePathname(filename, path=path,
                                                         mustExist=FALSE)
    verbose && cat(verbose, "Output pathname: ", pathname)

    if (!force && isFile(pathname)) {
      dfOut <- newInstance(clazz, filename=pathname)
      if (nbrOfUnits != nbrOfUnits(dfOut)) {
        throw("The number of units in existing output file does not match the number of units in the output file: ", nbrOfUnits, " != ", nbrOfUnits(dfOut))
      }
      verbose && cat(verbose, "Skipping already existing output file.")
      res[[kk]] <- dfOut
      verbose && exit(verbose)
      next
    }

    res[[kk]] %<-% {
      verbose && print(verbose, df)

      # Preallocate vector
      M <- rep(NA_real_, times=nbrOfUnits)

      verbose && enter(verbose, "Reading and smoothing input data")
      for (cc in seq_along(targetList)) {
        target <- targetList[[cc]]
        chromosome <- target$chromosome
        chrTag <- sprintf("Chr%02d", chromosome)
        verbose && enter(verbose, sprintf("Chromosome %d ('%s') of %d",
                                                 cc, chrTag, length(targetList)))

        verbose && cat(verbose, "Extracting raw CNs:")
        rawCNs <- extractRawCopyNumbers(df, chromosome=chromosome,
                                                    verbose=less(verbose, 10))
        verbose && print(verbose, rawCNs)
        verbose && summary(verbose, rawCNs)

        verbose && cat(verbose, "Smoothing CNs:")
        verbose && cat(verbose, "Target positions:")
        verbose && str(verbose, target$xOut)

        smoothCNs <- smoothRawCopyNumbers(this, rawCNs=rawCNs,
                                          target=target, verbose=verbose)
        rawCNs <- NULL; # Not needed anymore
        verbose && print(verbose, smoothCNs)
        verbose && summary(verbose, smoothCNs)

        M[target$units] <- getSignals(smoothCNs)

        target <- smoothCNs <- NULL # Not needed anymore

        verbose && exit(verbose)
      } # for (cc ...)
      verbose && exit(verbose)

      verbose && enter(verbose, "Smoothed CNs across all chromosomes:")
      verbose && str(verbose, M)
      verbose && summary(verbose, M)
      verbose && printf(verbose, "Missing values: %d (%.1f%%) out of %d\n",
                     sum(is.na(M)), 100*sum(is.na(M))/nbrOfUnits, nbrOfUnits)
      verbose && exit(verbose)

      verbose && enter(verbose, "Storing smoothed data")
      verbose && cat(verbose, "Pathname: ", pathname)

      footer <- list(
        sourceDataFile=list(
          fullname=getFullName(df),
          platform=getPlatform(df),
          chipType=getChipType(df),
          checksum=getChecksum(df)
        ),
        parameters=parameters
      )

      # Write to a temporary file
      pathnameT <- pushTemporaryFile(pathname, verbose=verbose)

      ## Alocate empty file
      dfOut <- clazz$allocate(filename=pathnameT, nbrOfRows=nbrOfUnits,
                              platform=parameters$targetUgp$platform,
                              chipType=parameters$targetUgp$chipType,
                              footer=footer, verbose=less(verbose, 50))

      ## Save signals
      dfOut[,1L] <- M

      # Not needed anymore
      M <- footer <- NULL

      # Renaming temporary file
      pathname <- popTemporaryFile(pathnameT, verbose=verbose)

      verbose && exit(verbose) # Storing

      # Sanity check
      dfOut <- newInstance(clazz, filename=pathname)
      if (nbrOfUnits != nbrOfUnits(dfOut)) {
        throw("The number of units in existing output file does not match the number of units in the output file: ", nbrOfUnits, " != ", nbrOfUnits(dfOut))
      }

      verbose && print(verbose, dfOut)

      dfOut
    } ## %<-%

    verbose && exit(verbose)
  } ## for (kk ...)

  ## Resolve futures
  res <- resolve(res)

  verbose && exit(verbose)

  dsOut <- getOutputDataSet(this, onMissing="error")

  invisible(dsOut)
}) # process()



setMethodS3("extractRawCopyNumbers", "BamDataFile", function(this, chromosome, ..., verbose=FALSE) {
  use("GenomicRanges")

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
# 2013-11-24
# o BUG FIX: smoothRawCopyNumbers() for TotalCnBinnedCounting assumed
#   that the 'matrixStats' package was loaded.
# 2012-10-11
# o Created from TotalCnBinnedSmoothing.R.
############################################################################

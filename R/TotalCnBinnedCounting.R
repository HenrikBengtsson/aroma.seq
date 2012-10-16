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
# @author
#*/########################################################################### 
setConstructorS3("TotalCnBinnedCounting", function(..., .reqSetClass="BamDataSet") {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  extend(TotalCnSmoothing(..., .reqSetClass=.reqSetClass), "TotalCnBinnedCounting");
})


setMethodS3("getExpectedOutputFullnames", "TotalCnBinnedCounting", function(this, ...) {
  names <- NextMethod("getExpectedOutputFullnames", this, ...);
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
  getPlatform(getFile(this, 1));
})

setMethodS3("getChipType", "BamDataSet", function(this, ...) {
  getChipType(getFile(this, 1));
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

###########################################################################/**
# @RdocClass TopHat2Alignment
#
# @title "The TopHat2Alignment class"
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
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
#  \item{indexSet}{An @see "Bowtie2IndexSet".}
#  \item{geneModelFile}{Gene model (transcriptome) GTF/GFF3 file.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux and OSX [1].
# }
#
# @author "TT"
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/}
# }
#*/###########################################################################
setConstructorS3("TopHat2Alignment", function(..., indexSet=NULL, geneModelFile=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indexSet':
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }

  # Argument 'geneModelFile':
  if (!is.null(geneModelFile)) {
    geneModelFile <- Arguments$getReadablePathname(geneModelFile);
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "TopHat2Alignment",
    geneModelFile = geneModelFile
  )
})



setMethodS3("getParameters", "TopHat2Alignment", function(this, ...) {
  params <- NextMethod("getParameters");
  params <- c(params, getOptionalArguments(this, ...));
  params;
}, protected=TRUE)


setMethodS3("getRootPath", "TopHat2Alignment", function(this, ...) {
  "tophat2Data";
}, protected=TRUE)


setMethodS3("getPath", "TopHat2Alignment", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this)

  # Platform
  ds <- getInputDataSet(this);
  platform <- "Generic";

  # The full path
  path <- filePath(rootPath, fullname, platform);

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


setMethodS3("getSampleNames", "TopHat2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  sampleNames <- sub("_(1|R1)$", "", getFullNames(ds));
  sampleNames;
})


setMethodS3("getExpectedOutputPaths", "TopHat2Alignment", function(this, ...) {
  # Find all available output directories
  path <- getPath(this);
  sampleNames <- getSampleNames(this);
  paths <- file.path(path, sampleNames);
  paths;
}, protected=TRUE)


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, ...) {
  ## Find all possible existing output data files
  bams <- BamDataSet$byPath(path=getPath(this), pattern="accepted_hits.bam$", recursive=TRUE);

  # Special case
  if (length(bams) == 0L) {
    bam <- BamDataFile(NA_character_, mustExist=FALSE);
    bams <- newInstance(bams, list(bam));
  }

  # Get the sample names in the found output set
  sampleNamesExpected <- getSampleNames(this);
  sampleNames <- getPathnames(bams);
  sampleNames <- basename(dirname(sampleNames));
  idxs <- match(sampleNamesExpected, sampleNames);
  bams <- extract(bams, idxs);

  # Sanity check
  stopifnot(length(bams) == length(sampleNamesExpected));

  bams;
})


setMethodS3("isDone", "TopHat2Alignment", function(this, ...) {
  bams <- getOutputDataSet(this);
  all(sapply(bams, FUN=isFile));
})


setMethodS3("process", "TopHat2Alignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Test for non-compatible bowtie2 and tophat2 versions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verT <- attr(findTopHat2(), "version");
  verB <- attr(findBowtie2(), "version");
  bad <- (verT == "2.0.3" && verB == "2.1.0");
  if (bad) {
    throw(sprintf("Detected incompatible software installations. TopHat2 v%s is known to not work with Bowtie2 v%s.", verT, verB))
  }
  verT <- verB <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "TopHat2 alignment");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify samples to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force) {
    todo <- seq_along(ds);
  } else {
    bams <- getOutputDataSet(this, verbose=less(verbose, 1));
    todo <- which(!sapply(bams, FUN=isFile));
  }
  verbose && cat(verbose, "Number of samples to process: ", length(todo));

  # Already done?
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already processed.");
    verbose && print(verbose, bams);
    verbose && exit(verbose);
    return(bams);
  }

  # The subset to be processed
  ds <- extract(ds, todo);
  verbose && cat(verbose, "Input data set (to be processed):");
  verbose && print(verbose, ds);

  isPaired <- isPaired(ds);
  verbose && cat(verbose, "Paired-end analysis: ", isPaired);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);

  outPath <- getPath(this);
  verbose && cat(verbose, "Output directory: ", outPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup arguments for TopHat
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gmf <- this$geneModelFile;
  if (!is.null(gmf)) optionsVec <- c("G"=gmf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(dfR1, isPaired=FALSE, indexSet, optionsVec=NULL, outPath, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'dfR1':
    dfR1 <- Arguments$getInstanceOf(dfR1, "BamDataFile");

    # Argument 'isPaired':
    isPaired <- Arguments$getLogical(isPaired);

    # Argument 'indexSet':
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'outPath':
    outPath <- Arguments$getWritablePath(outPath);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    sampleName <- sub("_(1|R1)$", "", getFullName(dfR1));
    verbose && enter(verbose, "Sample name ", sQuote(sampleName));

    verbose && cat(verbose, "R1 FASTQ file:");
    verbose && print(verbose, dfR1);

    args <- list(
      bowtieRefIndexPrefix=getIndexPrefix(indexSet),
      reads1=getPathname(dfR1),
      reads2=NULL,
      outPath=file.path(outPath, sampleName),
      optionsVec=optionsVec
    );

    if (isPaired) {
      dfR2 <- getMateFile(dfR1);
      verbose && cat(verbose, "R2 FASTQ file:");
      verbose && print(verbose, dfR2);
      args$reads2 <- getPathname(dfR2);
    }

    verbose && cat(verbose, "Arguments passed to TopHat:");
    verbose && str(verbose, args);

    res <- do.call(tophat2, args=args, verbose=less(verbose, 1));

    verbose && exit(verbose);

    invisible(list(res=res));
  }, isPaired=isPaired, indexSet=is, optionsVec=optionsVec, outPath=getPath(this), skip=skip, verbose=verbose) # dsApply()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bams <- getOutputDataSet(this, verbose=less(verbose, 1));
  verbose && print(verbose, bams);

  # Sanity check
  stopifnot(all(sapply(bams, FUN=isFile)));

  verbose && exit(verbose);

  bams;
})


############################################################################
# HISTORY:
# 2013-11-01 [HB]
# o SPEEDUP: Parallized process() for TopHat2Alignment.
# o Now process() for TopHat2Alignment skips already processed samples.
# o Now process() for TopHat2Alignment should also work for single-end
#   reads as well as paired-end reads.
# 2013-10-31 [HB]
# o Now utilizing tophat2().
# o CLEANUP: Dropped non-used argument 'outputDir' from constructor.
# o CLEANUP: Dropping stray cut'n'paste code.
# 2013-10-30 [HB]
# o ROBUSTNESS: Now process() for TopHat2Alignment checks for known
#   incompatible versions of bowtie2 and tophat2.
# 2013-08-12 [TT]
# o Created from Bowtie2Alignment.R.
############################################################################

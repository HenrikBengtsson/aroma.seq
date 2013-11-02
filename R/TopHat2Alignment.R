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
#   This method is available on Linux, OSX, and Windows [1].
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
  if (length(bams) > 0L) {
    # Keep the subset that corresponds to the input data set
    sampleNames <- sapply(bams, FUN=getPathname);
    sampleNames <- sapply(bams, FUN=dirname);
    sampleNames <- sapply(bams, FUN=basename);
    idxs <- match(getSampleNames(this), sampleNames);
    bams <- extract(bams, idxs);
  }
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
  args <- list(
    bowtieRefIndexPrefix=getIndexPrefix(is),
    reads1=NA_character_,
    reads2=NULL,
    outPath=NA_character_,
    optionsVec=character(0L)
  );
  gmf <- this$geneModelFile;
  if (!is.null(gmf)) args$optionsVec <- c("G"=gmf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process sample by sample...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (ii in seq_along(ds)) {
    dfR1 <- getFile(ds, ii);
    sampleName <- sub("_(1|R1)$", "", getFullName(dfR1));
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, sampleName, length(ds)));

    args$reads1 <- getPathname(dfR1);
    if (isPaired) {
      dfR2 <- getMateFile(dfR1);
      args$reads2 <- getPathname(dfR2);
    }
    args$outPath <- file.path(outPath, sampleName);
    res <- do.call(tophat2, args=args);

    verbose && exit(verbose);
  } # for (ii ...)


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

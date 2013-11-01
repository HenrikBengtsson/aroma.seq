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
#  \item{geneModelFile}{Gene model (transcriptome) gtf/gff file.}
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
    geneModelFile <- Arguments$getReadablePath(geneModelFile);
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "TopHat2Alignment");
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


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, ...) {
  ## Find all existing output data files
  res <- BamDataSet$byPath(path=getPath(this), pattern="accepted_hits.bam$", recursive=TRUE)

  ## TODO: Assert completeness
  res;
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
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "TopHat2 alignment");
  verbose && cat(verbose, "Paired-end analysis: ", isPaired(ds));

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);

  sampleNames <- sub("_1$", "", getFullNames(ds));
  verbose && cat(verbose, "Sample names:");
  verbose && print(verbose, sampleNames);

  outPath <- getPath(this);
  verbose && cat(verbose, "Output directory: ", outPath);

  # Setup arguments for tophat()
  args <- list(
    bowtieRefIndexPrefix=getIndexPrefix(is),
    reads1=NA,
    reads2=NA,
    outPath=NA,
    optionsVec=c("G"=this$geneModelFile)
  );

  if (isPaired(ds)) {
    filePairs <- getFilePairs(ds)
    fnsR1 <- sapply(filePairs[,1], FUN=getPathname)
    fnsR2 <- sapply(filePairs[,2], FUN=getPathname)
    for (ii in seq_along(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, getName(df), length(ds)));

      args$reads1 <- fnsR1[ii];
      args$reads2 <- fnsR2[ii];
      args$outPath <- file.path(outPath, sampleNames[ii]);
      res <- do.call(tophat2, args=args);

      verbose && exit(verbose);
    } # for (ii ...)
  } else {
    throw("Not yet implemented!");
    fnsR1 <- sapply(filePairs[,1], FUN=getPathname)  # <= BUG: Won't work! /HB 2013-11-01
    for (ii in seq_along(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, getName(df), length(ds)));

      args$reads1 <- fnsR1[ii];
      args$outPath <- file.path(outPath, sampleNames[ii]);
      res <- do.call(tophat2, args=args);

      verbose && exit(verbose);
    } # for (ii ...)
  }

  res <- getOutputDataSet(this, verbose=less(verbose, 1));
  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
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

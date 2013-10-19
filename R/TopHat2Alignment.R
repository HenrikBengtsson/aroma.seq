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
#  \item{outputDir}{A placeholder for the output path (overwritten by the code at present).}
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

setConstructorS3("TopHat2Alignment", function(..., indexSet=NULL, outputDir=NULL, geneModelFile=NULL) {
  # Validate arguments
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }
  if (!is.null(outputDir)) {
    outputDir <- Arguments$getWritablePath(outputDir)
    ## [ Convention:  This should be a single path, below which there can be multiple per-sample dirs ]
  }

  if (!is.null(geneModelFile)) {
    geneModelFile <- Arguments$getReadablePath(geneModelFile)
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
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asTopHatParameters <- function(rg, ...) {
    if (isEmpty(rg)) {
      return(NULL);
    }

    # Validate
##    if (!hasID(rg)) {
##      throw("... requires that the SAM read group has an ID.");
##    }

    rgArgs <- asString(rg, fmtstr="%s:%s");
    rgArgs <- rgArgs[regexpr("^ID:", rgArgs) == -1];

    # Don't forget to put within quotation marks
    rgArgs <- sprintf("\"%s\"", rgArgs);

    rgArgs <- as.list(rgArgs);
    names(rgArgs) <- rep("rg", times=length(rgArgs));

    rgArgs <- c(list("rg-id"=asSamList(rg)$ID), rgArgs);

    rgArgs;
  } # asTopHatParameters()

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
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);
  indexPrefix <- getIndexPrefix(is);

  rgSet <- this$.rgSet;
  if (!is.null(rgSet)) {
    verbose && cat(verbose, "Assigning SAM read group:");
    verbose && print(verbose, rgSet);
    validate(rgSet);
  }

  this$outputDir <- getPath(this)
  outputDirs <- file.path(this$outputDir, sub("_1$", "", getFullNames(ds)))

  if (isPaired(ds)) {
    filePairs <- getFilePairs(ds)
    fnsR1 <- sapply(filePairs[,1], function(x) {getPathname(x)})
    fnsR2 <- sapply(filePairs[,2], function(x) {getPathname(x)})
    for (i in seq_along(ds))
    {
      res <- do.call(what=tophat, args=list(bowtieRefIndexPrefix=getIndexPrefix(is),
                                            reads1=fnsR1[i],
                                            reads2=fnsR2[i],
                                            outDir=outputDirs[i],
                                            optionsVec=c("G"=this$geneModelFile),
                                            command="tophat2"))
      ## DEBUG
      cat(i, "\n")
    }
  } else {
    fnsR1 <- sapply(ds, function(x) {getPathname(x)})
    for (i in seq_along(ds))
    {
      res <- do.call(what=tophat, args=list(bowtieRefIndexPrefix=getIndexPrefix(is),
                                            reads1=fnsR1[i],
                                            outDir=outputDirs[i],
                                            optionsVec=c("G"=this$geneModelFile),
                                            command="tophat2"))
    }
  }

  params <- getParameters(this);
  verbose && cat(verbose, "Additional tophat2 arguments:");
  verbose && str(verbose, params);
  if (is.list(ds)) {
    verbose && cat(verbose, "Number of input data file pairs: ", length(ds[[1]]));
  } else {
    verbose && cat(verbose, "Number of input data files: ", length(ds));
  }

  res <- getOutputDataSet(this, verbose=less(verbose, 1));
  verbose && exit(verbose);
  invisible(res);
})


############################################################################
# HISTORY:
# 2013-08-12
# o Created from Bowtie2Alignment.R
############################################################################

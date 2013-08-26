###########################################################################/**
# @RdocClass TopHatAlignment
#
# @title "The TopHatAlignment class"
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
#  \item{outputDirs}{A placeholder for a vector of output paths (overwritten by the code at present).}
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
# \author{Taku Tokuyasu}
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/}
# }
#*/###########################################################################

setConstructorS3("TopHatAlignment", function(..., indexSet=NULL, outputDirs=NULL) {
  # Validate arguments
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }
  if (!is.null(outputDirs)) {
    outputDirs <- Arguments$getWritablePath(outputDirs)
    ## [ NB this is meant to be a vector in the case of multiple samples ]
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "TopHatAlignment");
})


setMethodS3("getParameters", "TopHatAlignment", function(this, ...) {
  params <- NextMethod("getParameters");
  params <- c(params, getOptionalArguments(this, ...));
  params;
}, protected=TRUE)


setMethodS3("getPath", "TopHatAlignment", function(this, create=TRUE, ...) {
  # [Create] Return the (sub-)directory tree for the data set
  
  # The full path
  # path <- filePath(rootPath, fullname, platform);
  path <- "tophat_out"  # Hard code the TopHat default for now
  
  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }
  
  # Verify that it is not the same as the input path
  if (is.list(ds)) {
    inPath <- getPath(ds[[1]])
  } else {
    inPath <- getPath(ds);
  }
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }
  
  path;
}, protected=TRUE)

setMethodS3("getOutputDataSet", "TopHatAlignment", function(this, ...) {
  ## Find all existing output data files
  res <- 
    sapply(this$outputDirs, function(dir) {
      BamDataSet$byPath(path=dir)
    }, simplify=FALSE)
  
  ## TODO: Assert completeness
  res;
})



setMethodS3("process", "TopHatAlignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
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
  
  verbose && enter(verbose, "TopHat alignment");
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

  # The following assumes ds is a list (i.e. paired-end data)
  fullNames <- sub("_1.fastq", "", basename(getPathnames(ds$read1)))
  outputDirs <- paste0("tophat_", fullNames)
  for (i in seq_along(ds$read1))
  {
    res <- do.call(what=tophat, args=list(getIndexPrefix(is),
                                          reads1=getPathnames(ds$read1)[i],
                                          reads2=getPathnames(ds$read2)[i],
                                          optionsVec=c("o"=outputDirs[i])))
  }

  this$outputDirs <- outputDirs  ## [ TAT - This is probably not aroma style]

  params <- getParameters(this);
  verbose && cat(verbose, "Additional tophat2 arguments:");
  verbose && str(verbose, params);
  if (is.list(ds)) { 
    verbose && cat(verbose, "Number of input data file pairs: ", length(ds[[1]]));
  } else {
    verbose && cat(verbose, "Number of input data files: ", length(ds));
  }
  # outPath <- getPath(this); # [ 201308 Punt on this for now; need to resolve getName, etc. for paired end datasets ]
  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2013-08-12
# o Created from Bowtie2Alignment.R
############################################################################

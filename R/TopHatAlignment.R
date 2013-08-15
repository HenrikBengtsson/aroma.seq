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

setConstructorS3("TopHatAlignment", function(..., indexSet=NULL) {
  # Validate arguments
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
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

  res <- do.call(what=tophat, args=list(getIndexPrefix(is),
                                        reads1=sapply(ds$read1, getPathname),
                                        reads2=sapply(ds$read2, getPathname)))
  
  params <- getParameters(this);
  verbose && cat(verbose, "Additional tophat2 arguments:");
  verbose && str(verbose, params);

  verbose && cat(verbose, "Number of files: ", length(ds));
  outPath <- getPath(this);
  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2013-08-12
# o Created from Bowtie2Alignment.R
############################################################################

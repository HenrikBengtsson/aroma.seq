###########################################################################/**
# @RdocClass QDNAseqEstimation
#
# @title "The QDNAseqEstimation class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "BamDataSet".}
#  \item{binWidth}{A positive @numeric specifying the bin width (in units
#    of kbp).
#    Alternatively, a @see "Biobase::AnnotatedDataFrame" specifying the bins.}
#  \item{log}{If @TRUE, the copy numbers are calculated on the log2 scale.}
#  \item{mappability, blacklist, residual, bases}{Post-filter arguments.}
#  \item{...}{Additional arguments passed to the QDNAseq method.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Benchmarking}{
#   TBA.
# }
#
# \references{
#  [1] TBA.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("QDNAseqEstimation", function(dataSet=NULL, binWidth=NULL, log=TRUE, mappability=50, blacklist=0, residual=2, bases=0, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");

    # Argument 'binWidth':
    if (inherits(binWidth, "AnnotatedDataFrame")) {
    } else {
      binWidth <- Arguments$getInteger(binWidth, range=c(0.1, 10e3));
    }
  } # if (!is.null(dataSet))

  extend(AromaSeqTransform(dataSet, binWidth=binWidth, log=log, mappability=mappability, blacklist=blacklist, residual=residual, bases=bases, ...), "QDNAseqEstimation")
})


setMethodS3("getBinWidth", "QDNAseqEstimation", function(this, ...) {
  # To please 'R CMD check'
  chromosome <- NULL; rm(list="chromosome");

  params <- getParameters(this);
  binWidth <- params$binWidth;

  if (inherits(binWidth, "AnnotatedDataFrame")) {
    bins <- binWidth;
    # AD HOC: Get 'binWidth' from loci
    data <- bins@data[,c("chromosome", "start")];
    chr <- data$chromosome[1L];
    data <- subset(data, chromosome == chr);
    binWidth <- median(diff(data$start), na.rm=TRUE);
    binWidth <- binWidth / 1e3;
    binWidth <- Arguments$getInteger(binWidth);
  }

  binWidth;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "QDNAseqEstimation", function(this, collapse=NULL, ...) {
  binWidth <- getBinWidth(this);
  binWidthTag <- sprintf("%gkb", binWidth);
  binWidthTag;
}, protected=TRUE)


setMethodS3("getRootPath", "QDNAseqEstimation", function(this, ...) {
  "qdnaseqData";
}, protected=TRUE)


setMethodS3("getOutputDataSet", "QDNAseqEstimation", function(this, ...) {
  ## Find all existing output data files
  path <- getPath(this);
  res <- RdsFileSet$byPath(path);
  ## Order according to input data set
  ds <- getInputDataSet(this);
  fullnames <- getFullNames(ds);
  res <- extract(res, fullnames, onMissing="NA");
  res;
}, protected=TRUE) # getOutputDataSet() for QDNAseqEstimation



setMethodS3("process", "QDNAseqEstimation", function(this, ..., force=FALSE, verbose=FALSE) {
  R.utils::use("QDNAseq (>= 1.2.4)");
  getBinAnnotations <- NULL; # To please 'R CMD check'.


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Disable 'QDNAseq' messages?
  if (!as.logical(verbose)) {
    oopts <- options("QDNAseq::verbose"=FALSE);
    on.exit(options(oopts), add=TRUE);
  }


  verbose && enter(verbose, "QDNAseq copy-number estimation");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  todo <- findFilesTodo(this, force=force, verbose=less(verbose, 1));
  # Already done?
  if (length(todo) == 0L) {
    verbose && cat(verbose, "Already done. Skipping.");
    res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
    verbose && exit(verbose);
    return(invisible(res));
  }
  verbose && cat(verbose, "Number of files to process: ", length(todo));

  params <- getParameters(this);
  verbose && cat(verbose, "Parameters:");
  verbose && str(verbose, params);

  binWidth <- params$binWidth;
  if (!inherits(binWidth, "AnnotatedDataFrame")) {
    verbose && enter(verbose, "QDNAseq/Retrieve QDNAseq bin annotation");
    bins <- getBinAnnotations(binWidth);
    verbose && print(verbose, bins);
    verbose && exit(verbose);
    binWidth <- bins;
  }
  params$binWidth <- binWidth;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds[todo], FUN=function(df, params, ..., force=FALSE, path, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq, QDNAseq (>= 0.5.6)");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

    # Argument 'params'
    params <- Arguments$getInstanceOf(params, "list");

    # Argument 'force':
    force <- Arguments$getLogical(force);

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "QDNAseq on one sample");

    # Arguments to doQDNAseq()
    args <- c(params, path=path, force=force);
    verbose && cat(verbose, "Arguments to doQDNAseq():");
    verbose && str(verbose, args);

    args <- c(list(df), args, verbose=verbose);
    cn <- do.call(doQDNAseq, args);
    verbose && print(verbose, cn);
    verbose && exit(verbose);

    invisible(list(cn=cn));
  }, params=params, path=getPath(this), force=force, verbose=verbose) # dsApply()

  # At this point, all files should have been processed
  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2013-11-16
# o Created.
############################################################################

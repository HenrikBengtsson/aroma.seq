###########################################################################/**
# @RdocDefault doQDNAseq
# @alias doQDNAseq.BamDataFile
# @alias doQDNAseq.BamDataSet
# @alias doQDNAseq.FastqDataSet
#
# @title "Quantitative inference of copy number aberrations with DNA isolated from fresh or formalin-fixed tissues by shallow whole-genome sequencing (QDNAseq)"
#
# \description{
#  @get "title" based on [1].
#  The algorithm is processed in bounded memory, meaning virtually
#  any number of samples can be analyzed on also very limited computer
#  systems.
# }
#
# \usage{
#   @usage doQDNAseq,FastqDataSet
#   @usage doQDNAseq,BamDataSet
#   @usage doQDNAseq,BamDataFile
# }
#
# \arguments{
#  \item{dataSet, df}{A @see "FastqDataSet" or a @see "BamDataSet" (or a @see "BamDataFile".}
#  \item{binWidth}{A positive @numeric specifying the bin width (in units of kbp).
#    Alternatively, a @see "Biobase::AnnotatedDataFrame" specifying the bins.}
#  \item{reference}{A @see "FastaReferenceFile" or a @see "BwaIndexSet" specifying the genome reference to align the FASTQ reads to.}
#  \item{log}{If @TRUE, the copy numbers are calculated on the log2 scale.}
#  \item{mappability, blacklist, residual, bases}{Post-filter arguments.}
#  \item{...}{Ignored, or passed to \code{doQDNAseq()}.}
#  \item{force}{If @TRUE, cached results are ignored.}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "R.filesets::GenericDataFileSet" containing
#   @see "QDNAseq::QDNAseqReadCounts" objects.
# }
#
# \references{
#  [1] TBA.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("doQDNAseq", "BamDataFile", function(df, binWidth, log=TRUE, mappability=50, blacklist=0, residual=2, bases=0, ..., path=".", force=FALSE, verbose=FALSE) {
  R.utils::use("QDNAseq (>= 0.5.8)");
  getBinAnnotations <- binReadCounts <- correctBins <- normalizeBins <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'binWidth':
  bins <- NULL;
  if (inherits(binWidth, "AnnotatedDataFrame")) {
    bins <- binWidth;
  } else {
    binWidth <- Arguments$getInteger(binWidth, range=c(0.1, 10e3));
  }

  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNAseq");
  verbose && print(verbose, df);

  # Output pathname
  filename <- sprintf("%s.RData", getFullName(df));
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  verbose && cat(verbose, "Output pathname: ", pathname);

  # Already done?
  isDone <- (!force && isFile(pathname));
  if (isDone) {
    verbose && cat(verbose, "Already processed. Skipping.");
    df <- GenericDataFile(pathname);
    verbose && exit(verbose);
    return(df);
  }


  # Disable 'QDNAseq' messages?
  if (!as.logical(verbose)) {
    oopts <- options("QDNAseq::verbose"=FALSE);
    on.exit(options(oopts), add=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(bins)) {
    verbose && enter(verbose, "QDNAseq/Retrieve QDNAseq bin annotation");
    bins <- getBinAnnotations(binWidth);
    verbose && print(verbose, bins);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNAseq/Reading and binning data");
  pathnameBAM <- getPathname(df);
  data <- binReadCounts(bins, bamfiles=pathnameBAM, cache=TRUE, force=force);
  verbose && print(verbose, data);
  bins <- pathnameBAM <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Correcting bin counts for GC content and mappability");
  dataC <- correctBins(data);
  verbose && print(verbose, dataC);
  data <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Normalization bin copy numbers");
  dataN <- normalizeBins(dataC, logTransform=log);
  verbose && print(verbose, dataN);
  dataC <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Saving");
  verbose && cat(verbose, "Output pathname: ", pathname);
  saveObject(dataN, file=pathname);
  verbose && exit(verbose);

  df <- GenericDataFile(pathname);

  verbose && exit(verbose);

  df;
}) # doQDNAseq()





setMethodS3("doQDNAseq", "BamDataSet", function(dataSet, ..., force=FALSE, verbose=FALSE) {
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


  verbose && enter(verbose, "QDNAseq");

  verbose && enter(verbose, "QDNAseq/copy number estimation");
  qe <- QDNAseqEstimation(dataSet, ...);
  verbose && print(verbose, qe);
  cns <- process(qe, force=force, verbose=verbose);
  verbose && print(verbose, cns);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cns;
}) # doQDNAseq()


setMethodS3("doQDNAseq", "FastqDataSet", function(dataSet, binWidth, reference, ..., verbose=FALSE) {
  R.utils::use("QDNAseq (>= 0.5.8)");
  getBinAnnotations <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'binWidth':
  bins <- NULL;
  if (inherits(binWidth, "AnnotatedDataFrame")) {
    bins <- binWidth;
  } else {
    binWidth <- Arguments$getInteger(binWidth, range=c(0.1, 10e3));
  }

  # Argument 'reference':
  if (inherits(reference, "FastaReferenceFile")) {
  } else if (inherits(reference, "BwaIndexSet")) {
  } else {
    throw("Argument 'reference' should either be of class 'FastaReferenceFile' or 'BwaIndexSet': ", class(reference)[1L]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNAseq");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNAseq/Check requirements");

  verbose && enter(verbose, "QDNAseq/Check requirements/BWA");
  stopifnot(isCapableOf(aroma.seq, "bwa"));
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Check requirements/Picard");
  stopifnot(isCapableOf(aroma.seq, "picard"));
  verbose && exit(verbose);

  verbose && exit(verbose);


  # Disable 'QDNAseq' messages?
  if (!as.logical(verbose)) {
    oopts <- options("QDNAseq::verbose"=FALSE);
    on.exit(options(oopts), add=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(bins)) {
    verbose && enter(verbose, "QDNAseq/Retrieve QDNAseq bin annotation");
    bins <- getBinAnnotations(binWidth);
    verbose && print(verbose, bins);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BWA alignment
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNAseq/BWA alignment");

  bs <- doBWA(dataSet, reference=reference, n=2, q=40, verbose=verbose);
  verbose && print(verbose, bs);

  # Not needed anymore
  reference <- NULL;

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove duplicated reads using Picard
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNAseq/Remove duplicated reads");
  dr <- PicardDuplicateRemoval(bs);
  verbose && print(verbose, dr);

  bsU <- process(dr, verbose=verbose);
  verbose && print(verbose, bsU);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # QDNAseq copy number estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNAseq/copy number estimation");
  cns <- doQDNAseq(bsU, binWidth=bins, ..., verbose=verbose);
  verbose && print(verbose, cns);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cns;
}) # doQDNAseq()


setMethodS3("doQDNAseq", "default", function(...) {
  throw("The \"default\" method is still not implemented. Please see help('doQDNAseq').");
})


############################################################################
# HISTORY:
# 2013-11-16
# o CLEANUP: Now doQDNAseq() for BamDataSet utilized QDNAseqEstimation.
# 2013-10-31
# o Updated capitalization to reflect the updated 'QDNAseq' package name.
# 2013-08-31
# o Now doQDNAseq for BamDataSet utilizes dsApply().
# 2013-08-31
# o BUG FIX: doQDNAseq() for BamDataSet would give an error when it
#   tried to collect and return the result file set.
# o BUG FIX: doQDNAseq() for BamDataSet would give an error if data set
#   was already processed and verbose output was enabled.
# 2013-08-22
# o CLEANUP: Now doQDNAseq() utilizes doBWA().
# 2013-07-29
# o Now doQDNAseq() for BamDataSet saves processed data to QDNAseqData/.
# 2013-07-11
# o SPEEDUP: Now doQDNAseq() only retrieves the bin annotation data
#   onces per call/data set.
# o Added Rdoc comments for doQDNAseq().
# o Added doQDNAseq() for FastqDataSet, which leverages ditto for
#   BamDataSet.
# 2013-07-03
# o Added to verbose statements.
# 2013-07-02
# o Created.
############################################################################

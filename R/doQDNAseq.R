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
#   Returns a @see "QDNAseq::QDNAseqReadCounts" object.
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
setMethodS3("doQDNAseq", "BamDataFile", function(df, binWidth, log=TRUE, mappability=50, blacklist=0, residual=2, bases=0, ..., force=FALSE, verbose=FALSE) {
  require("Biobase") || throw("Package not loaded: Biobase"); # combine()
  pkgName <- "QDNAseq";
  require(pkgName, character.only=TRUE) || throw("Package not loaded: QDNAseq");
  stopifnot(packageVersion("QDNAseq") >= "0.5.6");
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

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNAseq");
  verbose && print(verbose, df);


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
  pathname <- getPathname(df);
  data <- binReadCounts(bins, bamfiles=pathname, cache=TRUE, force=force);
  verbose && print(verbose, data);
  bins <- NULL;
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

  verbose && exit(verbose);

  dataN;
}) # doQDNAseq()





setMethodS3("doQDNAseq", "BamDataSet", function(dataSet, binWidth, ..., force=FALSE, verbose=FALSE) {
  pkgName <- "QDNAseq";
  require(pkgName, character.only=TRUE) || throw("Package not loaded: QDNAseq");
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

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNAseq");
  verbose && print(verbose, dataSet);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "QDNAseqData";
  path <- file.path(rootPath, getFullName(dataSet), "Generic");
  path <- Arguments$getWritablePath(path);


  if (is.null(bins)) {
    verbose && enter(verbose, "QDNAseq/Retrieve QDNAseq bin annotation");
    bins <- getBinAnnotations(binWidth);
    verbose && print(verbose, bins);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process each sample
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(dataSet, FUN=function(df, bins, ..., path, force=FALSE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

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

    verbose && enter(verbose, "QDNAseq of one sample");

    filename <- sprintf("%s.RData", getFullName(df));
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
    done <- (!force && isFile(pathname));
    if (done) {
      verbose && cat(verbose, "Already processed. Skipping.");
    } else {
      dataN <- doQDNAseq(df, binWidth=bins, ..., verbose=less(verbose,1));
      verbose && print(verbose, dataN);
      saveObject(dataN, file=pathname);
    } # if (done)

    verbose && exit(verbose);

    invisible(list(pathname=pathname));
  }, bins=bins, path=path, force=force, verbose=verbose) # dsApply()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Collect results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- GenericDataFileSet$byPath(path, pattern="[.]RData$");
  ds <- extract(ds, indexOf(ds, getNames(dataSet)));
  verbose && print(verbose, ds);

  verbose && exit(verbose);

  ds;
}) # doQDNAseq()


setMethodS3("doQDNAseq", "FastqDataSet", function(dataSet, binWidth, reference, ..., verbose=FALSE) {
  pkgName <- "QDNAseq";
  require(pkgName, character.only=TRUE) || throw("Package not loaded: QDNAseq");
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

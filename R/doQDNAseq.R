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
#  \item{...}{Additional arguments passed to @see "QDNAseq::applyFilters",
#    @see "QDNAseq::correctBins" and @see "QDNAseq::normalizeBins".}
#  \item{force}{If @TRUE, cached results are ignored.}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "R.filesets::RdsFileSet" containing
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
setMethodS3("doQDNAseq", "BamDataFile", function(df, binWidth, residual=TRUE, blacklist=TRUE, mappability=NA, bases=NA, filterAllosomes=TRUE, ..., path=".", force=FALSE, verbose=FALSE) {
  R.utils::use("QDNAseq (>= 1.2.4)");
  # To please 'R CMD check'
  getBinAnnotations <- binReadCounts <- applyFilters <- correctBins <- normalizeBins <- NULL;
  rm(list=c("getBinAnnotations", "binReadCounts", "applyFilters", "correctBins", "normalizeBins"), inherits=FALSE);

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
  filename <- sprintf("%s.rds", getFullName(df));
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  verbose && cat(verbose, "Output pathname: ", pathname);

  # Already done?
  isDone <- (!force && isFile(pathname));
  if (isDone) {
    verbose && cat(verbose, "Already processed. Skipping.");
    df <- RdsFile(pathname);
    verbose && exit(verbose);
    return(df);
  }


  # Disable 'QDNAseq' messages?
  if (!as.logical(verbose)) {
    oopts <- options("QDNAseq::verbose"=FALSE);
    on.exit(options(oopts), add=TRUE);
  }


  args <- list(...);
  keys <- names(args);

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
  data <- binReadCounts(bins, bamfiles=pathnameBAM, cache=FALSE, force=TRUE);
  verbose && print(verbose, data);
  bins <- pathnameBAM <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Filtering out bins based on a priori filters");
  keysF <- setdiff(names(formals(applyFilters)), c("object", "force", "..."))
  argsF <- args[intersect(keys, keysF)]
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, argsF);
  dataF <- do.call(applyFilters, args=c(list(data), argsF));
  verbose && print(verbose, dataF);
  data <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Correcting bin counts for GC content and mappability");
  keysC <- setdiff(names(formals(correctBins)), c("object", "force", "..."))
  argsC <- args[intersect(keys, keysC)];
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, argsC);
  dataC <- do.call(correctBins, args=c(list(dataF), argsC));
  verbose && print(verbose, dataC);
  dataF <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Normalization bin copy numbers");
  keysN <- setdiff(names(formals(normalizeBins)), c("object", "force", "..."))
  argsN <- args[intersect(keys, keysN)];
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, argsN);
  dataN <- do.call(normalizeBins, args=c(list(dataC), argsN));
  verbose && print(verbose, dataN);
  dataC <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Saving");
  verbose && cat(verbose, "Output pathname: ", pathname);
  saveRDS(dataN, file=pathname);
  verbose && exit(verbose);

  df <- RdsFile(pathname);

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
  stop_if_not(isCapableOf(aroma.seq, "bwa"));
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNAseq/Check requirements/Picard");
  stop_if_not(isCapableOf(aroma.seq, "picard"));
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

  # Daud Sie (VUMC) confirm that they are using 'bwa aln -n 2 -q 40'
  # [private email from HB on 2013-11-19]
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

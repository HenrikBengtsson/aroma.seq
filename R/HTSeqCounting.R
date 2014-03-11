###########################################################################/**
# @RdocClass HTSeqCounting
#
# @title "The HTSeqCounting class"
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
#  \item{dataSet}{A @see "BamDataSet".}
#  \item{transcripts}{A @see "GtfDataFile" specifying a gene model.}
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   ...
# }
#
# @author "HB"
#
# \references{
#  [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
#      data with Python}, EMBL, Jan 2014.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr
# }
#*/###########################################################################
setConstructorS3("HTSeqCounting", function(dataSet=NULL, transcripts=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");

    # Argument 'transcripts':
    transcripts <- Arguments$getInstanceOf(transcripts, "GtfDataFile");
  } # if (!is.null(dataSet))


  # Arguments '...':
  args <- list(...);

  extend(AromaSeqTransform(dataSet, ...), "HTSeqCounting",
    .transcripts = transcripts
  );
})


setMethodS3("getRootPath", "HTSeqCounting", function(this, ...) {
  "htseqCountData";
}, protected=TRUE)


setMethodS3("getParameters", "HTSeqCounting", function(this, ...) {
  params <- NextMethod("getAsteriskTags");
  params$transcripts <- this$.transcripts;
  params;
}, protected=TRUE)


setMethodS3("getSampleNames", "HTSeqCounting", function(this, ...) {
  ds <- getInputDataSet(this);
  getFullNames(ds, ...);
}, protected=TRUE)


setMethodS3("process", "HTSeqCounting", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert external software versions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stopifnot(isCapableOf(aroma.seq, "htseq"));


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

  verbose && enter(verbose, "HTSeq counting");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify groups to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force) {
    todo <- seq_along(ds);
  } else {
    counts <- getOutputDataSet(this, onMissing="NA", verbose=less(verbose, 1));
    todo <- which(!sapply(counts, FUN=isFile));
  }
  verbose && cat(verbose, "Number of samples to process: ", length(todo));

  # Already done?
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already processed.");
    verbose && print(verbose, counts);
    verbose && exit(verbose);
    return(counts);
  }

  outPath <- getPath(this);
  verbose && cat(verbose, "Output directory: ", outPath);

  params <- getParameters(this);
  transcripts <- params$transcripts;
  verbose && cat(verbose, "Using transcripts:");
  verbose && print(verbose, transcripts);
  # Workaround for *gzipped* GTF files (not supported by HTSeq binaries)
  if (isGzipped(transcripts)) {
    verbose && enter(verbose, "Temporary uncompressing file");
    pathnameZ <- getPathname(transcripts)
    pathname <- gunzip(pathnameZ, temporary=TRUE, remove=FALSE)
    on.exit(file.remove(pathname), add=TRUE);
    transcripts <- newInstance(transcripts, pathname);
    verbose && cat(verbose, "Using (temporary) transcripts:");
    verbose && print(verbose, transcripts);
    verbose && exit(verbose);
  }
  # Sanity check
  stopifnot(!isGzipped(transcripts));

  verbose && cat(verbose, "Number of samples: ", length(ds));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating compatibility between GTF and BAM data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Validating compatibility between GTF and BAM data set");
  # Get one BAM file.
  bam <- ds[[1L]];
  bamNames <- getTargetNames(bam);
  verbose && printf(verbose, "Target/chromosome names in BAM: %s [%d]\n", hpaste(bamNames), length(bamNames));

  # Slightly faster to parse unzipped file.
  gtf <- transcripts;
  gtfNames <- getSeqNames(gtf, unique=TRUE);
  gtf <- params$transcripts;
  verbose && printf(verbose, "Sequence/chromosome names in GTF: %s [%d]\n", hpaste(gtfNames), length(gtfNames));

  inBoth <- intersect(bamNames, gtfNames);
  if (length(inBoth) == 0L) {
    msg <- sprintf("Incompatible GTF file: None of the %d sequence/chromosome names in the GTF (%s) are found in the BAM file (%s): [%s] not in [%s]", length(gtfNames), getPathname(gtf), getPathname(bam), hpaste(gtfNames), hpaste(bamNames));
    verbose && cat(verbose, msg);
    throw(msg);
  }

  inBAMnotGTF <- setdiff(bamNames, gtfNames);
  if (length(inBAMnotGTF) > 0L) {
    msg <- sprintf("Possibly an incompatible GTF file: Found %d target/chromosome names in the BAM file (%s) that are not in the GTF (%s): [%s] not in [%s]", length(inBAMnotGTF), getPathname(bam), getPathname(gtf), hpaste(inBAMnotGTF), hpaste(gtfNames));
    verbose && cat(verbose, msg);
    # FIXME: Should we use throw() here instead?  Set some thresholding? /HB 2014-03-10
    warning(msg);
  }

  inGTFnotBAM <- setdiff(gtfNames, bamNames);
  if (length(inGTFnotBAM) > 0L) {
    msg <- sprintf("Not counting all sequences/chromosomes in GTF file: Found %d sequence/chromosome names in the GTF (%s) that are not in the BAM file (%s): [%s] not in [%s]", length(inGTFnotBAM), getPathname(gtf), getPathname(bam), hpaste(inGTFnotBAM), hpaste(gtfNames));
    verbose && cat(verbose, msg);
    warning(msg);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, IDXS=todo, DROP=TRUE, FUN=function(df, transcripts=NULL, outPath, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq, Rsamtools");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

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

    verbose && enter(verbose, "htseq-count");

    verbose && cat(verbose, "BAM file:");
    verbose && print(verbose, df);
    pathnameBAM <- getPathname(df);

    pathnameGTF <- getPathname(transcripts);
    verbose && cat(verbose, "GTF pathname: ", pathnameGTF);

    filenameD <- sprintf("%s.count", getFullName(df));
    pathnameD <- Arguments$getWritablePathname(filenameD, path=outPath, mustNotExist=FALSE);
    # Nothing to do?
    if (skip && isFile(pathnameD)) {
      verbose && cat(verbose, "Already processed. Skipping.");
      verbose && exit(verbose);
      return(GenericDataFile(pathnameD));
    }

    # Final sample-specific output directory
    args <- list(
      pathnameS=pathnameBAM,
      gff=pathnameGTF,
      pathnameD=pathnameD
    );

    # Is BAM file sorted?
    isSorted <- isSorted(df);
    if (isSorted) {
      # ...then assume it is sorted by position (aroma.seq policy)
      # FIXME: Check how BAM file is sorted
      args$orderedBy <- "position";
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BEGIN: ATOMIC OUTPUT
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##    # Write to temporary output directory
##    args$outPath <- sprintf("%s.tmp", args$outPath);
##    verbose && cat(verbose, "Temporary output directory: ", args$outPath);

    verbose && cat(verbose, "Arguments passed to htseqCount():");
    verbose && str(verbose, args);
    args$verbose <- less(verbose, 1);
    res <- do.call(htseqCount, args=args);

    verbose && cat(verbose, "Results:");
    verbose && print(verbose, res);

    # Was there a non-zero exit status?
    status <- attr(res, "status");
    if (!is.null(status)) {
      verbose && cat(verbose, "Status: ", status);
    }

##   # Rename from temporary to final directory
##    file.rename(args$outPath, outPathS);
##    verbose && cat(verbose, "Final output directory: ", outPathS);
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # END: ATOMIC OUTPUT
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    verbose && exit(verbose);

    invisible(list(res=res));
  }, transcripts=transcripts, outPath=getPath(this), skip=skip, verbose=verbose) # dsApply()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  counts <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
  verbose && print(verbose, counts);

  # Sanity check
  stopifnot(all(sapply(counts, FUN=isFile)));

  verbose && exit(verbose);

  counts;
})


############################################################################
# HISTORY:
# 2014-03-10 [HB]
# o Now process() passes orderedBy="position" to htseqCount().
# o ROBUSTNESS: Now process() for HTSeqCounting validates that the GTF
#   file has chromsome names in a format that matched the BAM files.
#   If not, an error is thrown.
# 2014-01-23 [HB]
# o Created.
############################################################################

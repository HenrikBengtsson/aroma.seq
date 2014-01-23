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
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
#  \item{transcripts}{A @see "GtfDataFile" specifying a gene model.}
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
setConstructorS3("HTSeqCounting", function(..., transcripts=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'transcripts':
  transcripts <- Arguments$getInstanceOf(transcripts, "GtfDataFile");

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(...), "HTSeqCounting",
    transcripts = transcripts
  )
})


setMethodS3("getRootPath", "HTSeqCounting", function(this, ...) {
  "htseqCountData";
}, protected=TRUE)


setMethodS3("getParameters", "HTSeqCounting", function(this, ...) {
  params <- NextMethod("getAsteriskTags");
  params$transcripts <- this$transcripts;
  params;
}, protected=TRUE)


setMethodS3("getSampleNames", "HTSeqCounting", function(this, ...) {
  ds <- getInputDataSet(this);
  getFullNames(ds, ...);
}, protected=TRUE)

setMethodS3("getExpectedOutputPaths", "HTSeqCounting", function(this, ...) {
  # Find all available output directories
  path <- getPath(this);
  sampleNames <- getSampleNames(this);
  paths <- file.path(path, sampleNames);
  paths;
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

  transcripts <- this$transcripts;
  verbose && cat(verbose, "Using transcripts:");
  verbose && print(verbose, transcripts);
##   # Workaround for *gzipped* GTF files (not supported by HTSeq binaries)
##   if (isGzipped(transcripts)) {
##     verbose && enter(verbose, "Temporary uncompressing file");
##     pathnameZ <- getPathname(transcripts)
##     pathname <- gunzip(pathnameZ, temporary=TRUE, remove=FALSE)
##     on.exit(file.remove(pathname), add=TRUE);
##     transcripts <- newInstance(transcripts, pathname);
##     verbose && cat(verbose, "Using (temporary) transcripts:");
##     verbose && print(verbose, transcripts);
##     verbose && exit(verbose);
##   }
  # Sanity check
  stopifnot(!isGzipped(transcripts));

  verbose && cat(verbose, "Number of samples: ", length(ds));


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

    pathnameBAM <- getPathname(df);
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

    pathnameGTF <- getPathname(transcripts);
    verbose && cat(verbose, "GTF pathname: ", pathnameGTF);

    filenameD <- sprintf("%s.count", getFullName(df));
    pathnameD <- Arguments$getWritablePathname(filenameD, path=outPath, mustNotExist=FALSE);
    # Nothing to do?
    if (skip && isFile(pathname)) {
      verbose && cat(verbose, "Already processed. Skipping.");
      verbose && exit(verbose);
      return(GenericDataFile(pathnameD));
    }

    # Final sample-specific output directory
    outPathS <- file.path(outPath, sampleName);
    args <- list(
      pathnameS=pathnameBAM,
      gff=pathnameGTF
      pathnameD=pathnameD
    );


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BEGIN: ATOMIC OUTPUT
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##    # Write to temporary output directory
##    args$outPath <- sprintf("%s.tmp", args$outPath);
##    verbose && cat(verbose, "Temporary output directory: ", args$outPath);

    # Run htseq-count
    verbose && cat(verbose, "Arguments passed to htseq-count:");
    verbose && str(verbose, args);
    args$verbose <- less(verbose, 1);
    res <- do.call(htseqCount, args=args);

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
# 2014-01-23 [HB]
# o Created.
############################################################################

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
#  \item{transcripts}{Gene model (transcriptome) GTF/GFF3 file.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux and OSX [1].
# }
#
# @author "TT"
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/} \cr
#  [2] Trapnell et al. \emph{Differential gene and transcript expression
#      analysis of RNA-seq experiments with TopHat and Cufflinks}.
#      Nat Protoc, 2012.\cr
# }
#*/###########################################################################
setConstructorS3("TopHat2Alignment", function(..., indexSet=NULL, transcripts=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'indexSet':
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }

  # Argument 'transcripts':
  if (!is.null(transcripts)) {
    transcripts <- Arguments$getInstanceOf(transcripts, "GenericDataFile");
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "TopHat2Alignment",
    transcripts = transcripts
  )
})


setMethodS3("getRootPath", "TopHat2Alignment", function(this, ...) {
  "tophat2Data";
}, protected=TRUE)


setMethodS3("getSampleNames", "TopHat2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  sampleNames <- sub("_(1|R1)$", "", getFullNames(ds));
  sampleNames;
})

setMethodS3("getExpectedOutputPaths", "TopHat2Alignment", function(this, ...) {
  # Find all available output directories
  path <- getPath(this);
  sampleNames <- getSampleNames(this);
  paths <- file.path(path, sampleNames);
  paths;
}, protected=TRUE)


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
  # Argument 'force':
  force <- Arguments$getLogical(force);

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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify samples to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force) {
    todo <- seq_along(ds);
  } else {
    bams <- getOutputDataSet(this, onMissing="NA", verbose=less(verbose, 1));
    todo <- which(!sapply(bams, FUN=isFile));
  }
  verbose && cat(verbose, "Number of samples to process: ", length(todo));

  # Already done?
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already processed.");
    verbose && print(verbose, bams);
    verbose && exit(verbose);
    return(bams);
  }

  # The subset to be processed
  ds <- extract(ds, todo);
  verbose && cat(verbose, "Input data set (to be processed):");
  verbose && print(verbose, ds);

  isPaired <- isPaired(ds);
  verbose && cat(verbose, "Paired-end analysis: ", isPaired);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);

  outPath <- getPath(this);
  verbose && cat(verbose, "Output directory: ", outPath);

  transcripts <- this$transcripts;
  verbose && cat(verbose, "Using transcripts:");
  verbose && print(verbose, transcripts);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(dfR1, isPaired=FALSE, indexSet, transcripts=NULL, outPath, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'dfR1':
    dfR1 <- Arguments$getInstanceOf(dfR1, "FastqDataFile");

    # Argument 'isPaired':
    isPaired <- Arguments$getLogical(isPaired);

    # Argument 'indexSet':
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");

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

    sampleName <- sub("_(1|R1)$", "", getFullName(dfR1));
    verbose && enter(verbose, "Sample name ", sQuote(sampleName));

    gtf <- NULL;
    if (!is.null(transcripts)) gtf <- getPathname(transcripts);

    verbose && cat(verbose, "R1 FASTQ file:");
    verbose && print(verbose, dfR1);

    # Final sample-specific output directory
    outPathS <- file.path(outPath, sampleName);
    args <- list(
      bowtieRefIndexPrefix=getIndexPrefix(indexSet),
      reads1=getPathname(dfR1),
      reads2=NULL,
      outPath=outPathS,
      gtf=gtf
    );

    if (isPaired) {
      dfR2 <- getMateFile(dfR1);
      verbose && cat(verbose, "R2 FASTQ file:");
      verbose && print(verbose, dfR2);
      args$reads2 <- getPathname(dfR2);
    }

    verbose && cat(verbose, "Arguments passed to TopHat:");
    verbose && str(verbose, args);

    args$verbose <- less(verbose, 1);

    # BEGIN: Atomic output
    # Write to temporary output directory
    args$outPathS <- sprintf("%s.tmp", args$outPathS);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Align reads using TopHat2
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    res <- do.call(tophat2, args=args);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (b) Generates BAM index file (assuming the BAM file is sorted)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pathnameBAM <- file.path(args$outPathS, "accepted_hits.bam");
    pathnameBAI <- indexBam(pathnameBAM);

    # Rename from temporary to final directory
    file.rename(args$outPathS, outPathS);
    # END: Atomic output

    verbose && exit(verbose);

    invisible(list(res=res));
  }, isPaired=isPaired, indexSet=is, transcripts=transcripts, outPath=getPath(this), skip=skip, verbose=verbose) # dsApply()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bams <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
  verbose && print(verbose, bams);

  # Sanity check
  stopifnot(all(sapply(bams, FUN=isFile)));

  verbose && exit(verbose);

  bams;
})


############################################################################
# HISTORY:
# 2013-11-21 [HB]
# o Now process() for TopHat2Alignment generates an index file for
#   accepted_hits.bam.  This will also assert the assumption that
#   TopHat2 outputs sorted BAM files, because if not indexing will
#   fail and an error will be generated.
# 2013-11-16 [HB]
# o CLEANUP: Dropped several methods now taken care of by super class.
# 2013-11-01 [HB]
# o SPEEDUP: Parallized process() for TopHat2Alignment.
# o Now process() for TopHat2Alignment skips already processed samples.
# o Now process() for TopHat2Alignment should also work for single-end
#   reads as well as paired-end reads.
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

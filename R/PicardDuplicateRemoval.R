###########################################################################/**
# @RdocClass PicardDuplicateRemoval
#
# @title "The PicardDuplicateRemoval class"
#
# \description{
#  @classhierarchy
#
#  This method \emph{flags} reads that are aligned to more than one locus,
#  which is done using Picard's 'MarkDuplicates' method [1].
#
#  Note that it is assumed that the input BAM files are already sorted,
#  which also means that it can be assumed that the output BAM files
#  are sorted.  As with all other methods that outputs BAM files,
#  this methods index all outputted BAM files.
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "BamDataSet".}
#  \item{...}{Additional arguments passed to @see "AromaSeqTransform".}
#  \item{ASSUME_SORTED, VALIDATION_STRINGENCY}{
#    Additional arguments passed to Picard MarkDuplicates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Benchmarking}{
#  As a very rough guideline, a 1.0GB BAM file takes
#  about 10-15 minutes to process using this method.
# }
#
# \references{
#  [1] Picard, \url{http://picard.sourceforge.net/}\cr
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("PicardDuplicateRemoval", function(dataSet=NULL, ..., ASSUME_SORTED=TRUE, VALIDATION_STRINGENCY="SILENT") {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");
  } # if (!is.null(dataSet))

  extend(AromaSeqTransform(dataSet, ASSUME_SORTED=ASSUME_SORTED, VALIDATION_STRINGENCY=VALIDATION_STRINGENCY, ...), "PicardDuplicateRemoval");
})


setMethodS3("getAsteriskTags", "PicardDuplicateRemoval", function(this, collapse=NULL, ...) {
  "-dups";
}, protected=TRUE)



setMethodS3("getRootPath", "PicardDuplicateRemoval", function(this, ...) {
  # Use same root path as input data set
  ds <- getInputDataSet(this);
  path <- getPath(ds);
  path <- getParent(path, depth=2L);
  # Sanity check
  stopifnot(regexpr("Data$", path) != -1L);
  path;
}, protected=TRUE)



setMethodS3("process", "PicardDuplicateRemoval", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Picard removal of duplicated reads");

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  if (force) {
    todo <- seq_along(ds);
  } else {
    todo <- findFilesTodo(this, verbose=less(verbose, 1));
    # Already done?
    if (length(todo) == 0L) {
      verbose && cat(verbose, "Already done. Skipping.");
      res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
      verbose && exit(verbose);
      return(invisible(res));
    }
  }

  nbrOfFiles <- length(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  params <- getParameters(this);
  verbose && cat(verbose, "Additional Picard MarkDuplicates arguments:");
  verbose && str(verbose, params);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds[todo], FUN=function(df, params, path, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "Picard MarkDuplicates on one sample");

    pathname <- getPathname(df);
    verbose && cat(verbose, "Input BAM pathname: ", pathname);

    # Output BAM file
    fullname <- getFullName(df);
    filename <- sprintf("%s.bam", fullname);
    pathnameD <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "Output BAM pathname: ", pathnameD);
    pathnameDI <- gsub("[.]bam$", ".bai", pathnameD);
    verbose && cat(verbose, "Output BAM index pathname: ", pathnameDI);

    # Nothing to do?
    done <- (skip && isFile(pathnameD) && isFile(pathnameDI));
    if (done) {
      verbose && cat(verbose, "Already processed. Skipping");
    } else {
      verbose && print(verbose, df);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Filter via Picard
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!skip || !isFile(pathnameD)) {
        verbose && enter(verbose, "Calling Picard MarkDuplicates");

        pathnameM <- gsub("[.]bam$", ",picard,MarkDuplicates,metrics.log", pathnameD);
        args <- list(
          "MarkDuplicates",
          INPUT=pathname,
          OUTPUT=pathnameD,
          METRICS_FILE=pathnameM,
          REMOVE_DUPLICATES=TRUE,
          ASSUME_SORTED=TRUE,            # TODO: Assert! /HB 2012-10-02
          VALIDATION_STRINGENCY="SILENT"
        );
        verbose && cat(verbose, "Arguments:");
        verbose && str(verbose, df);

        # Assert no overwrite
        stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));
        stopifnot(getAbsolutePath(pathnameM) != getAbsolutePath(pathnameD));

        args$verbose <- less(verbose, 20);
        res <- do.call(systemPicard, args);
        status <- attr(res, "status"); if (is.null(status)) status <- 0L;
        verbose && cat(verbose, "Results:");
        verbose && str(verbose, res);
        verbose && cat(verbose, "Status:");
        verbose && str(verbose, status);

        verbose && exit(verbose);
      } # if (!isFile(pathnameD))

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Generic BAM index
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      bf <- BamDataFile(pathnameD);
      verbose && print(verbose, bf);

      if (!skip || !hasIndex(bf)) {
        verbose && enter(verbose, "Creating BAM index");
        bfi <- buildIndex(bf, skip=skip, overwrite=!skip, verbose=less(verbose, 10));
        verbose && exit(verbose);
      }
    } # if (done)

    verbose && exit(verbose);

    invisible(list(pathnameD=pathnameD, pathnameDI=pathnameDI));
  }, params=params, path=getPath(this), skip=!force, verbose=verbose) # dsApply()

  # At this point, all files should have been processed
  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));

  verbose && exit(verbose);

  res;
})



############################################################################
# HISTORY:
# 2013-11-16
# o CLEANUP: Dropped several methods now taken care of by super class.
# 2013-11-15
# o Added argument 'onMissing' to getOutputDataSet().
# 2013-09-03
# o Now process() for PicardDuplicateRemoval utilizes dsApply().
# 2012-11-26
# o BUG FIX: getOutputDataSet() would return a data set with "missing"
#   files, if not complete.  Now it only returns the existing files.
# 2012-10-02
# o Created.
############################################################################

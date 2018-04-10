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
#  \item{...}{Arguments passed to @see "AromaSeqTransform".}
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
setConstructorS3("PicardDuplicateRemoval", function(..., ASSUME_SORTED=TRUE, VALIDATION_STRINGENCY="SILENT") {
  extend(SamTransform(..., ASSUME_SORTED=ASSUME_SORTED, VALIDATION_STRINGENCY=VALIDATION_STRINGENCY), "PicardDuplicateRemoval");
})


setMethodS3("getAsteriskTags", "PicardDuplicateRemoval", function(this, collapse=NULL, ...) {
  "-dups";
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
  path <- getPath(this)
  skip <- !force

  res <- listenv()
  for (kk in seq_along(todo)) {
    ii <- todo[kk]
    df <- ds[[ii]]
    fullname <- getFullName(df)

    verbose && enter(verbose, sprintf("Picard MarkDuplicates on sample #%d ('%s') of %d", kk, fullname, length(todo)))

    pathname <- getPathname(df)
    verbose && cat(verbose, "Input BAM pathname: ", pathname)

    # Output BAM file
    filename <- sprintf("%s.bam", fullname)
    pathnameD <- Arguments$getWritablePathname(filename, path=path)
    verbose && cat(verbose, "Output BAM pathname: ", pathnameD)
    pathnameDI <- gsub("[.]bam$", ".bai", pathnameD)
    verbose && cat(verbose, "Output BAM index pathname: ", pathnameDI)

    # Nothing to do?
    done <- (skip && isFile(pathnameD) && isFile(pathnameDI))
    if (done) {
      verbose && cat(verbose, "Already processed. Skipping")
      res[[kk]] <- pathnameD
      verbose && exit(verbose)
      next
    }

    res[[kk]] %<-% {
      verbose && print(verbose, df)

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
        stop_if_not(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));
        stop_if_not(getAbsolutePath(pathnameM) != getAbsolutePath(pathnameD));

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
        buildIndex(bf, skip=skip, overwrite=!skip, verbose=less(verbose, 10));
        verbose && exit(verbose);
      }

      pathnameD <- Arguments$getReadablePathname(pathnameD)
      pathnameD
    } ## %<-%

    verbose && exit(verbose)
  } ## for (kk ...)

  ## Resolve all futures
  res <- resolve(res)

  # At this point, all files should have been processed
  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1))

  verbose && exit(verbose)

  res
})

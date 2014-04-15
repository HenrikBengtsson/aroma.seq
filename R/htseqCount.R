###########################################################################/**
# @RdocDefault htseqCount
#
# @title "Calls the htseq-count executable to count input reads on features"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameS}{An input BAM or SAM file containing aligned reads.}
#   \item{gff}{The gene feature file, in GFF/GTF format}
#   \item{orderedBy}{A @character string specifying how the input file has
#    been sorted, if at all.}
#   \item{sortByName}{A @character string specifying when the BAM/SAM file
#    should be sorted by name.}
#   \item{optionsVec}{A named @character @vector of options to htseq-count.}
#   \item{...}{(Not used)}
#   \item{pathnameD}{(optional) destination file to save htseq-count output.}
#   \item{command}{A @character string specifying the name of the executable.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns what @see "systemHTSeqCount" returns.
# }
#
# \section{Backward compatibility}{
#   \code{htseq-count} (< 0.6.0) requires (i) a SAM file as input
#   that (ii) is sorted by name.  However, this method will take of that
#   internally, iff needed.  That is, it will created a temporary SAM
#   file that is sorted by query name before passing it to \code{htseq-count}.
# }
#
# \references{
#  [1] S Anders, TP Pyl, W Huber,
#      HTSeq - A Python framework to work with high-throughput sequencing data.
#      bioRxiv 2014. doi: 10.1101/002824.\cr
# }
#
# @author "TT,HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("htseqCount", "default", function(pathnameS, gff, orderedBy=c("none", "position", "name"),
                                              sortByName=c("always", "auto"),
                                              optionsVec=c('-s'='no', '-a'='10'),  ## Anders et al default
                                              ...,
                                              pathnameD=NULL, command='htseq-count', verbose=FALSE) {
  ## ( Support a call like this: "htseq-count -s no -a 10 lib_sn.sam gff > countFile")
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  R.utils::use("Rsamtools")

  isSAM <- function(pathname, ...) {
    (regexpr(".*[.]sam$", pathname, ignore.case=TRUE) != -1L);
  } # isSAM()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameS'
  pathnameS <- Arguments$getReadablePathname(pathnameS);
  # Assert BAM or SAM file
  if (regexpr(".*[.](bam|sam)$", pathnameS, ignore.case=TRUE) == -1L) {
    throw("Not a BAM or SAM file: ", pathnameS);
  }

  # Argument 'gff'
  gff <- Arguments$getReadablePathname(gff);

  # Argument 'orderedBy':
  orderedBy <- match.arg(orderedBy);

  # Argument 'sortByName':
  sortByName <- match.arg(sortByName);

  # Count file
  if (is.null(pathnameD)) {
    pathnameD <- sub("[.](bam|sam)$", ".count", pathnameS, ignore.case=TRUE);
  }
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Running htseqCount()");

  # Log file
  pathnameL <- sprintf("%s.log", pathnameD);
  pathnameL <- Arguments$getWritablePathname(pathnameL, mustNotExist=TRUE);


  verbose && cat(verbose, "Input BAM/SAM file: ", sQuote(pathnameS));
  verbose && cat(verbose, "Genome transcript file: ", sQuote(gff));
  verbose && cat(verbose, "Output count file: ", sQuote(pathnameD));
  verbose && cat(verbose, "Output log file: ", sQuote(pathnameL));
  verbose && cat(verbose, "htseq-count executable: ", command);

  # Locates the htseq-count executable
  bin <- findHTSeq(command=command, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);
  ver <- attr(bin, "version");
  ver <- numeric_version(ver);
  verbose && cat(verbose, "Version: ", ver);
  if (is.na(ver)) ver <- numeric_version("0.0.0");


  # BACKWARD COMPATIBILITY for htseq-count (< 0.6.0)
  convertToSam <- "auto";
  if (convertToSam == "auto") {
    if (ver < "0.6.0") {
      verbose && cat(verbose, "Backward compatibility: Detected htseq-count (< 0.6.0): Will convert to SAM [EXPENSIVE WORKAROUND]");
      convertToSam <- "always";
    } else {
      convertToSam <- "never";
    }
  }
  convertToSam <- match.arg(convertToSam, c("always", "never"));

  # BACKWARD COMPATIBILITY for htseq-count (< 0.6.0)
  if (sortByName == "auto") {
    if (orderedBy == "name") {
      sortByName <- "never";
    } else if (ver < "0.6.0") {
      verbose && cat(verbose, "Backward compatibility: Detected htseq-count (< 0.6.0): Will sort BAM/SAM by name [EXPENSIVE WORKAROUND]");
      sortByName <- "always";
    } else {
      if (orderedBy == "none") {
        sortByName <- "always";
      } else {
        sortByName <- "never";
      }
    }
  }
  sortByName <- match.arg(sortByName, c("always", "never"));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sort by query name?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (sortByName == "always") {
    verbose && enter(verbose, "Making sure to provide with a BAM/SAM file sorted by query name [EXPENSIVE WORKAROUND]");

    # (a) Create/setup BAM input for sorting, iff SAM file is passed
    if (isSAM(pathnameS)) {
      verbose && enter(verbose, "Converting BAM to SAM");
      pathnameBAM <- sub("[.]sam$", ",tmp.bam", pathnameS, ignore.case=TRUE);
      pathnameBAM <- Arguments$getWritablePathname(pathnameBAM, mustNotExist=TRUE);
      on.exit({
        if (isFile(pathnameBAM)) file.remove(pathnameBAM);
      }, add=TRUE);
      # Convert SAM-to-BAM
      asSam(pathnameS, pathnameBAM);
      # Sanity check
      pathnameBAM <- Arguments$getReadablePathname(pathnameBAM);
      verbose && cat(verbose, "Temporary BAM file: ", pathnameBAM);
      # Update input pathname
      pathnameS <- pathnameBAM;
      verbose && exit(verbose);
    } else {
      pathnameBAM <- NULL;
    }

    # (b) Sort input BAM by name
    verbose && enter(verbose, "Sorting BAM by query name");
    pathnameBAMs <- sub("(|,tmp)[.]bam$", ",byName,tmp.bam", pathnameS);
    pathnameBAMs <- Arguments$getWritablePathname(pathnameBAMs, mustNotExist=TRUE);
    on.exit({
      if (isFile(pathnameBAMs)) file.remove(pathnameBAMs);
    }, add=TRUE);
    # Sort BAM
    sortBam(pathnameS, sub("[.]bam$", "", pathnameBAMs), byQname=TRUE)
    # Sanity check
    pathnameBAMs <- Arguments$getReadablePathname(pathnameBAMs);
    verbose && cat(verbose, "Sorted temporary BAM file: ", pathnameBAMs);
    # Remove temporary non-sorted BAM file
    if (isFile(pathnameBAM)) file.remove(pathnameBAM);
    # Update input pathname
    pathnameS <- pathnameBAMs;
    orderedBy <- "name";
    verbose && exit(verbose);

    verbose && exit(verbose);
  } else {
    pathnameBAMs <- NULL;
  } # if (sortByName == "always")

  # Sanity check
  pathnameS <- Arguments$getReadablePathname(pathnameS);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Convert BAM to SAM?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (convertToSam == "always" && !isSAM(pathnameS)) {
    verbose && enter(verbose, "Converting sorted BAM to SAM");

    pathnameSAMs <- sub("[.]bam$", ".sam", pathnameS);
    pathnameSAMs <- Arguments$getWritablePathname(pathnameSAMs);
    on.exit({
      if (isFile(pathnameSAMs)) file.remove(pathnameSAMs);
    }, add=TRUE);
    # Convert BAM-to-SAM
    samtoolsView(pathnameS, pathnameSAMs);
    # Sanity checks
    pathnameSAMs <- Arguments$getReadablePathname(pathnameSAMs);
    verbose && cat(verbose, "Sorted temporary SAM file: ", pathnameSAMs);
    # Update input pathname
    pathnameS <- pathnameBAMs;

    # Remove temporary sorted BAM file
    if (isFile(pathnameBAMs)) file.remove(pathnameBAMs);

    verbose && exit(verbose);
  }

  # Sanity check
  pathnameS <- Arguments$getReadablePathname(pathnameS);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call htseq-count
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling systemHTSeqCount()");
  opts <- NULL;
  if (ver >= "0.6.0") {
    # SAM or BAM format?
    if (isSAM(pathnameS)) {
      opts <- c(opts, "--format=sam");
    } else {
      opts <- c(opts, "--format=bam");
    }

    # Ordered by position or by name?
    if (orderedBy == "position") {
      opts <- c(opts, "--order=pos");
    } else if (orderedBy == "name") {
      opts <- c(opts, "--order=name");
    } else {
      throw("Invalid value on argument 'orderedBy': ", orderedBy);
    }
  } # if (ver >= "0.6.0")

  # Add user options
  opts <- c(opts, optionsVec);

  # All htseq-count arguments including pathnames
  args <- c(opts, shQuote(pathnameS), shQuote(gff));

  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);

  # htseq-count sends results to standard output, so make sure
  # to redirect to output file.

  # Write results and log to temporary files
  pathnameDT <- pushTemporaryFile(pathnameD, verbose=verbose);
  pathnameLT <- pushTemporaryFile(pathnameL, verbose=verbose);

  status <- systemHTSeqCount(..., args=args, stdout=pathnameDT, stderr=pathnameLT, command=command, verbose=less(verbose, 1));
  verbose && cat(verbose, "Exit status: ", status);

  # Renaming temporary log file
  pathnameL <- popTemporaryFile(pathnameLT, verbose=verbose);
  pathnameL <- Arguments$getReadablePathname(pathnameL);


  # Was there a non-zero exit status?
  if (is.numeric(status) && status != 0L) {
    verbose && enter(verbose, "Non-zero exit status: ", status);

    # Remove empty count file
    fi <- file.info(pathnameDT);
    if (fi$size == 0L) {
      file.remove(pathnameDT);
      verbose && cat(verbose, "Removed empty result file: ", pathnameDT);
      pathnameDT <- NULL;
    }

    # Parse log file
    log <- readLines(pathnameL);

    # Did an error occur?  Then throw an R error
    idx <- grep("^Error", log, value=FALSE, ignore.case=TRUE)[1L];
    if (!is.na(idx)) {
      msg <- log[seq(from=idx,to=length(log))];
      msg <- paste(msg, collapse="\n");
      msg <- sprintf("Error detected while running %s (v%s) on %s [exit code %d]: %s", command, ver, shQuote(pathnameS), status, msg);
      throw(msg);
    }

    # Was a warning generated?  The translate it into an R warning
    idx <- grep("^Warning", log, value=FALSE, ignore.case=TRUE)[1L];
    if (!is.na(idx)) {
      msg <- log[seq(from=idx,to=length(log))];
      msg <- paste(msg, collapse="\n");
      msg <- sprintf("Warning detected while running %s (v%s) on %s [exit code %d]: %s", command, ver, shQuote(pathnameS), status, msg);
      warning(msg);
    }
    verbose && exit(verbose);
  }

  # Renaming temporary count file
  if (!is.null(pathnameDT)) {
    pathnameD <- popTemporaryFile(pathnameDT, verbose=verbose);
    pathnameD <- Arguments$getReadablePathname(pathnameD);
  }

  res <- character(0L);
  attr(res, "status") <- status;
  attr(res, "countFile") <- pathnameD;
  attr(res, "logFile") <- pathnameL;

  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # htseqCount()


############################################################################
# HISTORY:
# 2014-03-11 [HB]
# o ROBUSTNESS: Added argument 'sortByName' to htseqCount(), which for
#   now defaults to "always", because although htseq-count (>= 0.6.0)
#   is supposed to handle when BAM files are sorted by position, it will
#   run out of memory for modestly large BAM files.
# o ROBUSTNESS: If the external software returns a non-zero exit code,
#   then htseqCount() removes any empty count file and scans the log file
#   for errors and warnings and translates them into ditto in R.
# o Now htseqCount() writes counts to one file and log messages to another.
# 2014-03-10 [HB]
# o ROBUSTNESS: Now htseqCount() uses shQuote() for all pathnames.
# 2014-03-09 [HB]
# o ROBUSTNESS: Now htseqCount() writes results atomically.
# o CLEANUP: Now htseqCount() cleans out temporary files ASAP.
# o SPEEDUP: Now supporting htseq-count (>= 0.6.0) options.
# o CLEANUP: Clean out obsolete code handling 'optionsVec'.
# o Added verbose output to htseqCount().
# o Temporary filenames are now using aroma-style tags.
# o DOCUMENTATION: Added Details section to help.
# 2013-12-16 [TT]
# o local asSam() added (needed for Rsamtools < 1.15.14).
# 2013-11-29 [TT]
# o Internal sort-by-name added.
# 2013-05-31 [TT]
# o Created.
############################################################################

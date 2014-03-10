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
#   \item{pathnameD}{(optional) destination file to save htseq-count output.}
#   \item{gff}{The gene feature file, in GFF/GTF format}
#   \item{orderedBy}{A @character string specifying how the input file has
#    been sorted, if at all.}
#   \item{optionsVec}{A named @character @vector of options to htseq-count.}
#   \item{...}{(Not used)}
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
setMethodS3("htseqCount", "default", function(pathnameS,
                                              pathnameD=NULL,
                                              gff,
					      orderedBy=c("none", "position", "name"),
                                              optionsVec=c('-s'='no', '-a'='10'),  ## Anders et al default
                                              ...,
                                              command='htseq-count',
                                              verbose=FALSE) {

  ## ( Support a call like this: "htseq-count -s no -a 10 lib_sn.sam gff > countFile")
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  R.utils::use("Rsamtools")

  # BACKWARD COMPATIBILITY: Add asSam(), iff missing.
  if (packageVersion("Rsamtools") < "1.15.14") {
    asSam <- function(file, destination, ...) {
      file <- Arguments$getReadablePathname(file)
      fileD <- sprintf("%s.sam", destination)
      fileD <- Arguments$getWritablePathname(fileD)
      samtoolsView(file, fileD)
      fileD
    } # asSam()
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameS'
  pathnameS <- Arguments$getReadablePathname(pathnameS);
  # Assert BAM or SAM file
  if (regexpr(".*[.](bam|sam)$", pathnameS, ignore.case=TRUE) == -1L) {
    throw("Not a BAM or SAM file: ", pathnameS);
  }

  # Argument 'pathnameD'
  if (is.null(pathnameD)) {
    pathnameD <- sub("[.](bam|sam)$", ".count", pathnameS, ignore.case=TRUE)
    if (pathnameD == pathnameS) { # should not happen
      pathnameD <- paste(pathnameS, ".count", sep="")
    }
  }
  stopifnot(pathnameD != pathnameS)
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=TRUE);

  # Argument 'gff'
  gff <- Arguments$getReadablePathname(gff);

  # Argument 'orderedBy':
  orderedBy <- match.arg(orderedBy);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Running htseqCount()");
  verbose && cat(verbose, "Input BAM/SAM file: ", sQuote(pathnameS));
  verbose && cat(verbose, "Genome transcript file: ", sQuote(gff));
  verbose && cat(verbose, "Output count file: ", sQuote(pathnameD));
  verbose && cat(verbose, "htseq-count executable: ", command);

  # Locates the htseq-count executable
  bin <- findHTSeq(command=command, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);
  ver <- attr(bin, "version");
  ver <- numeric_version(ver);
  verbose && cat(verbose, "Version: ", ver);
  if (is.na(ver)) ver <- numeric_version("0.0.0");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BACKWARD COMPATIBILITY for htseq-count (< 0.6.0)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ver < "0.6.0") {
    verbose && enter(verbose, "Making sure to provide htseq-count (< 0.6.0) with a SAM file sorted by query name [EXPENSIVE WORKAROUND]");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Create/setup BAM input for sorting, iff SAM file is passed
    isSAM <- (regexpr(".*[.]sam$", pathnameS, ignore.case=TRUE) != -1L);
    if (isSAM) {
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
    pathnameBAMs <- sub("[.]bam$", ",byName,tmp.bam", pathnameS)
    pathnameBAMs <- Arguments$getWritablePathname(pathnameBAMs, mustNotExist=TRUE)
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
    verbose && exit(verbose);

    # (c) Convert sorted BAM to SAM
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
    # Remove temporary sorted BAM file
    if (isFile(pathnameBAMs)) file.remove(pathnameBAMs);
    # Update input pathname
    pathnameS <- pathnameBAMs;
    verbose && exit(verbose);

    verbose && exit(verbose);
  } else {
    # Sort input BAM?
    if (orderedBy == "none") {
      verbose && enter(verbose, "Sorting non-sorted BAM by query name");
      pathnameBAMs <- sub("[.]bam$", ",byName,tmp.bam", pathnameS);
      pathnameBAMs <- Arguments$getWritablePathname(pathnameBAMs, mustNotExist=TRUE);
      on.exit({
        if (isFile(pathnameBAMs)) file.remove(pathnameBAMs);
      }, add=TRUE);
      # Sort BAM
      sortBam(pathnameS, sub("[.]bam$", "", pathnameBAMs), byQname=TRUE);
      # Sanity check
      pathnameBAMs <- Arguments$getReadablePathname(pathnameBAMs);
      verbose && cat(verbose, "Sorted temporary BAM file: ", pathnameBAMs);
      # Update input pathname
      pathnameS <- pathnameBAMs;
      orderedBy <- "name";
      verbose && exit(verbose);
    }
  } # if (... ver < "0.6.0")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call htseq-count
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling systemHTSeqCount()");
  opts <- NULL;
  if (ver >= "0.6.0") {
    # SAM or BAM format?
    if (regexpr("[.]bam$", pathnameS, ignore.case=TRUE) != -1L) {
      opts <- c(opts, "--format=bam");
    } else {
      opts <- c(opts, "--format=sam");
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

  # Write to a temporary file
  pathnameDT <- pushTemporaryFile(pathnameD, verbose=verbose);
  args <- c(args, " > ", shQuote(pathnameDT));

  res <- do.call(what=systemHTSeqCount, args=list(command=command, args=args));

  # Renaming temporary file
  pathnameD <- popTemporaryFile(pathnameDT, verbose=verbose);

  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # htseqCount()


############################################################################
# HISTORY:
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

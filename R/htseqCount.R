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
# \details{
#   The \code{htseq-count} executable requires (i) a SAM file as input
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Make sure to provide htseq-count with a SAM file sorted by query name
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Create/setup BAM input for sorting, iff SAM file is passed
  isSAM <- (regexpr(".*[.]sam$", pathnameS, ignore.case=TRUE) != -1L);
  if (isSAM) {
    verbose && enter(verbose, "Converting BAM to SAM");
    inBam <- sub("[.]sam$", ",tmp.bam", pathnameS, ignore.case=TRUE);
    inBam <- Arguments$getWritablePathname(inBam, mustNotExist=TRUE);
    asSam(pathnameS, inBam);
    on.exit({
      file.remove(inBam)
    }, add=TRUE)
    # Sanity check
    inBam <- Arguments$getReadablePathname(inBam);
    verbose && exit(verbose);
  } else {
    inBam <- pathnameS;
  }

  # (b) Sort input BAM by name
  verbose && enter(verbose, "Sorting BAM by query name");
  inBamS <- sub("[.]bam$", ",byName,tmp.bam", inBam)
  inBamS <- Arguments$getWritablePathname(inBamS, mustNotExist=TRUE)
  sortBam(inBam, sub("[.]bam$", "", inBamS), byQname=TRUE)
  # Sanity check
  inBamS <- Arguments$getReadablePathname(inBamS)
  verbose && exit(verbose);

  # (c) Convert sorted BAM to SAM (needed by htseq-count)
  verbose && enter(verbose, "Converting sorted BAM to SAM");
  inSamS <- sub("[.]bam$", ".sam", inBamS)
  inSamS <- Arguments$getWritablePathname(inSamS)
  samtoolsView(inBamS, inSamS)
  on.exit({
    file.remove(inBamS, inSamS)
  }, add=TRUE)
  # Sanity checks
  inSamS <- Arguments$getReadablePathname(inSamS)
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call htseq-count
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling systemHTSeqCount()");
  htseqArgs <- c(inSamS, gff)
  htseqArgs <- c(htseqArgs, " > ", pathnameD)
  htseqOptions <- NULL
  if (!is.null(optionsVec)) {
    htseqOptions <- optionsVec
    # Deprecate this (dashes now required to be added by the user)
    # nms <- names(htseqOptions); names(htseqOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }
  args <- list(command=command, args=c(htseqOptions, htseqArgs));
  verbose && cat(verbose, "Arguments passed:");
  verbose && str(verbose, args);

  res <- do.call(what=systemHTSeqCount, args=args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2014-03-09 [HB]
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

###########################################################################/**
# @RdocFunction bowtie2
#
# @title "Calls the Bowtie2 executable to align input reads"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reads1}{(required) A @vector of FASTQ pathnames of reads.}
#   \item{reads2}{(optional; paired-end only) A @vector of FASTQ pathnames of mate reads.}
#   \item{indexPrefix}{Bowtie2 reference index prefix.}
#   \item{pathnameSAM}{Output SAM file.}
#   \item{...}{...}
#   \item{gzAllowed}{A @logical specifying whether gzipped FASTQ files are supported or not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns ...
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
bowtie2 <- function(reads1, reads2=NULL, indexPrefix, pathnameSAM, ..., gzAllowed=NA, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasCommas <- function(pathnames, ...) {
    (regexpr(",", pathnames, fixed=TRUE) != -1L);
  } # hasCommas()

  assertNoCommas <- function(pathnames, ...) {
    stopifnot(!any(hasCommas(pathnames)));
  } # assertNoCommas()

		       
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reads1'
  reads1 <- Arguments$getReadablePathnames(reads1, absolute=TRUE);
  assertNoDuplicated(reads1);

  # Argument 'reads2'
  isPaired <- (length(reads2) > 0L);
  if (isPaired) {
    stopifnot(length(reads2) == length(reads1));
    reads2 <- Arguments$getReadablePathnames(reads2, absolute=TRUE);
    assertNoDuplicated(reads2);
  } 

  # Argument 'indexPrefix':
  indexPrefix <- Arguments$getCharacter(indexPrefix);

  # Argument 'pathnameSAM':
  pathnameSAM <- Arguments$getWritablePathname(pathnameSAM);
  assertNoDuplicated(pathnameSAM);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Running bowtie2()");

  verbose && cat(verbose, "R1 FASTQ files:");
  verbose && print(verbose, sapply(reads1, FUN=getRelativePath));
  if (isPaired) {
     verbose && cat(verbose, "R2 FASTQ files:");
     verbose && print(verbose, sapply(reads2, FUN=getRelativePath));
  } 
  verbose && cat(verbose, "Paired alignment: ", isPaired);
  verbose && cat(verbose, "Reference index prefix: ", indexPrefix);
  outPath <- dirname(pathnameSAM);
  verbose && cat(verbose, "Output directory: ", outPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handle gzip'ed FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnameFQ <- c(reads1, reads2);
  assertNoDuplicated(pathnameFQ);
  isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
  if (isGzipped) {
    verbose && enter(verbose, "Detected gzipped FASTQ files");

    if (is.na(gzAllowed)) {
      gzAllowed <- queryBowtie2("support:fastq.gz");
    }

    if (!gzAllowed) {
      verbose && enter(verbose, "Temporarily decompressing gzipped FASTQ files");

      decompress <- getOption(aromaSettings, "devel/fastq.gz/decompress", TRUE);
      if (!decompress) {
        why <- attr(gzAllowed, "why");
        throw(sprintf("Cannot align reads in '%s': %s", getPathname(df), why));
      }

      # If not, temporarily decompress (=remove when done)
      pathnameFQ <- sapply(pathnameFQ, FUN=gunzip, temporary=TRUE, remove=FALSE);
      pathnameFQtmpA <- pathnameFQ;
      on.exit({
        # Make sure to remove temporary file
        lapply(pathnameFQtmpA, FUN=function(pathname) {
          if (isFile(pathname)) file.remove(pathname);
        });
      }, add=TRUE);

      # Sanity check
      isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
      stopifnot(!isGzipped);

      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  } # if (isGzipped)



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # WORKAROUND: Bowtie2() does not support commas in the FASTQ
  # pathname.  If so, use a temporary filename without commas.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasComma <- (regexpr(",", pathnameFQ, fixed=TRUE) != -1L);
  if (any(hasComma)) {
    pathnameFQ[hasComma] <- sapply(pathnameFQ[hasComma], FUN=function(pathname) {
      ext <- if (isGzipped(pathname)) ".fq.gz" else ".fq";
      pathnameT <- tempfile(fileext=ext);
      createLink(pathnameT, target=pathname);
      pathnameT;
    });

    # Remove temporary files
    on.exit({
      file.remove(pathnameFQ[hasComma]);
    }, add=TRUE);
  }

  # Split up in reads1 and reads2 again
  n <- length(reads1);
  reads1 <- pathnameFQ[seq_len(n)];
  if (isPaired) {
    reads2 <- pathnameFQ[-seq_len(n)];
    # Sanity check
    stopifnot(length(reads1) == length(reads2));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call Bowtie2 executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && exit(verbose, "Calling systemBowtie2()");

  indexPath <- Arguments$getReadablePath(getParent(indexPrefix));
  verbose && cat(verbose, "Index path: ", indexPath);
  args <- list("-x"=shQuote(indexPrefix))
  if (isPaired) {
    args[["-1"]] <- shQuote(reads1);
    args[["-2"]] <- shQuote(reads2);
  } else {
    args[["-U"]] <- shQuote(reads1);
  }
  args[["-S"]] <- shQuote(pathnameSAM);
  args <- c(args, list(...));
  verbose && cat(verbose, "Arguments:");
  verbose && print(verbose, args); 
  res <- systemBowtie2(args=args, verbose=verbose);
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status:");
  verbose && str(verbose, status); 

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # In case bowtie2 generates empty SAM files
  # /HB 2012-10-01 (still to be observed)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isFile(pathnameSAM)) {
    if (file.info(pathnameSAM)$size == 0L) {
      verbose && cat(verbose, "Removing empty SAM file falsely created by Bowtie2: ", pathnameSAM);
      file.remove(pathnameSAM);
    }
  }

  verbose && exit(verbose);

  res;
} # bowtie2()


############################################################################
# HISTORY:
# 2014-08-07 [HB]
# o Now using arguments 'reads1' and 'reads2', cf. tophat().
# o Added verbose output.
# 2014-03-10 [HB]
# o ROBUSTNESS: Now bowtie2() uses shQuote() for all pathnames.
# 2014-01-14 [HB]
# o ROBUSTNESS: Now bowtie2() tests for duplicated entries in 'reads1'
#   and 'reads2' and gives an informative errors message if detected.
# 2013-08-24
# o Now bowtie2() will do paired-end alignment if length(pathnameFQ) == 2.
# 2013-08-23
# o BUG FIX: Read Group options ('--rg' and '--rg-id') passed to 'bowtie2'
#   by the Bowtie2Aligment class missed the preceeding '--'.  Also, if
#   the Read Group ID was missing NULL was used - now it is set to 1.
# 2013-07-18
# o Now bowtie2() handles if there are commas in the pathname of
#   the FASTQ file by using a tempory file link without commas.  This
#   is needed because the bowtie2 executable does not support commas.
# 2013-06-27
# o Now bowtie2() temporarily decompresses gzipped FASTQ files in case
#   the installed bowtie2 does not support gzip files.
# 2012-09-27
# o Created.
############################################################################

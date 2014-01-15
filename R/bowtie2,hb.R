bowtie2_hb <- function(pathnameFQ, indexPrefix, pathnameSAM, ..., gzAllowed=NA, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameFQ':
  pathnameFQ <- Arguments$getReadablePathnames(pathnameFQ, length=1:2);
  assertNoDuplicated(pathnameFQ);

  # Argument 'indexPrefix':
  indexPrefix <- Arguments$getCharacter(indexPrefix);

  # Argument 'pathnameSAM':
  pathnameSAM <- Arguments$getWritablePathname(pathnameSAM);
  assertNoDuplicated(pathnameSAM);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handle gzip'ed FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
  if (isGzipped) {
    if (is.na(gzAllowed)) {
      gzAllowed <- queryBowtie2("support:fastq.gz");
    }

    if (!gzAllowed) {
      decompress <- getOption(aromaSettings, "devel/fastq.gz/decompress", TRUE);
      if (!decompress) {
        why <- attr(gzAllowed, "why");
        throw(sprintf("Cannot align reads in '%s': %s", getPathname(df), why));
      }

      # If not, temporarily decompress (=remove when done)
      pathnameFQ <- sapply(pathnameFQ, FUN=gunzip, temporary=TRUE, remove=FALSE);
      on.exit({
        # Make sure to remove temporary file
        lapply(pathnameFQ, FUN=function(pathname) {
          if (isFile(pathname)) file.remove(pathname);
        });
      }, add=TRUE);

      # Sanity check
      isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
      stopifnot(!isGzipped);
    }
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call external bowtie2 executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isPaired <- (length(pathnameFQ) == 2L);
  indexPath <- Arguments$getReadablePath(getParent(indexPrefix));
  if (isPaired) {
    # Sanity check
    stopifnot(pathnameFQ[1L] != pathnameFQ[2L]);
    res <- systemBowtie2(args=list("-x"=indexPrefix, "-1"=pathnameFQ[1L], "-2"=pathnameFQ[2L], "-S"=pathnameSAM, ...), verbose=verbose);
  } else {
    res <- systemBowtie2(args=list("-x"=indexPrefix, "-U"=pathnameFQ, "-S"=pathnameSAM, ...), verbose=verbose);
  }


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

  res;
} # bowtie2_hb()


############################################################################
# HISTORY:
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

###########################################################################/**
# @RdocDefault bowtie2Build
#
# @title "Creates index on reference genome using bowtie2-build"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameFAs}{A @character @vector of FASTA reference files.}
#   \item{bowtieRefIndexPrefix}{A @character string specifying the bowtie2
#     reference index to be built (partial pathname, minus the .*.bt2 suffix).}
#   \item{optionsVec}{(optional) A named @character @vector.}
#   \item{...}{...}
#   \item{command}{The name of the external executable.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \section{Support for compressed input files}{
#   If gzip'ed FASTA files are used, this function will temporarily decompress
#   before passing them to the bowtie2-build external software (which only
#   support non-compressed FASTA files).
# }
#
# \section{Known issues}{
#   The FASTA pathnames must not contain commas.
#   If detected, this method generates an informative error.
# }
#
# @author "HB,TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bowtie2Build", "default", function(pathnameFAs,
                                                bowtieRefIndexPrefix,
                                                optionsVec=NULL,
                                                ...,
                                                command="bowtie2-build",
                                                verbose=FALSE) {
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
  # Argument 'pathnameFAs':
  pathnameFAs <- Arguments$getReadablePathnames(pathnameFAs);
  assertNoDuplicated(pathnameFAs);
  assertNoCommas(pathnameFAs);

  # Argument 'bowtieRefIndexPrefix'
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix, length=c(1L,1L));
  bowtieRefIndexPath <- dirname(bowtieRefIndexPrefix);
  bowtieRefIndexPath <- Arguments$getWritablePath(bowtieRefIndexPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Workaround for gzip'ed FASTA files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isGzipped <- sapply(pathnameFAs, FUN=isGzipped);
  if (any(isGzipped)) {
    # Temporarily decompress gzip'ed FASTA files
    pathnameFAs[isGzipped] <- sapply(pathnameFAs[isGzipped], FUN=gunzip, remove=FALSE, temporary=TRUE);
    on.exit({
      file.remove(pathnameFAs[isGzipped]);
    }, add=TRUE);
  }


  # Sanity check
  assertNoCommas(pathnameFAs);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup call to bowtie2-build binary
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add dashes as appropriate to names of "bowtie2 options"
  opts <- optionsVec;
  if (length(opts) > 0L) {
    nms <- names(opts);
    names(opts) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="");
  }

  # Append FASTA reference files
  opts <- c(opts, paste(unname(pathnameFAs), collapse=","));

  # Append bowtie reference index prefix
  opts <- c(opts, bowtieRefIndexPrefix);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call bowtie2-build binary
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(command=command, args=opts);
  res <- do.call(what=systemBowtie2Build, args=args);

  res;
}) # bowtie2Build()

############################################################################
# HISTORY:
# 2014-01-14 [HB]
# o ROBUSTNESS: Now bowtie2Build() tests for duplicated entries in
#   'pathnameFAs' and gives an informative errors message if detected.
# 2013-11-01 [HB]
# o Now bowtie2Build() supports gzip'ed FASTQ.
# 2013-03-08 [TT]
# o Completely rewritten to follow tophat.R template.
# 2012-09-14 [TT]
# o First real draft of "level 2" code.
# 2012-07-18 [HB]
# o Created bowtie2-build() stub.
############################################################################

###########################################################################/**
# @RdocDefault bwaSamse
#
# @title "Generates BWA-backtrack single-end (SE) alignments via 'bwa samse'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameSAI}{The SAI file to be aligned.}
#   \item{pathnameFQ}{The FASTQ file to be aligned.}
#   \item{indexPrefix}{The pathname prefix to the BWA index files.}
#   \item{pathnameD}{The destination pathname.}
#   \item{...}{Additional arguments specifying BWA 'samse' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
#   bwaIndex(pathnameFA)
#
#   pathnameSAI <- "bwaData/LambdaVirusExample/Lambda_phage/reads_1.sai";
#   pathnameFQ <- "fastqData/LambdaVirusExample/Lambda_phage/reads_1.fq";
#   pathnameD <- "bwaData/LambdaVirusExample/Lambda_phage/reads_1.sam";
#   bwaSamse(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ,
#            pathnameFA=pathnameFA, pathnameD=pathnameD);
# }}
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaSamse", "default", function(pathnameSAI, pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameSAI':
  pathnameSAI <- Arguments$getReadablePathname(pathnameSAI);

  # Argument 'pathnameFQ':
  pathnameFQ <- Arguments$getReadablePathname(pathnameFQ);

  # Argument 'indexPrefix':
  dummy <- Arguments$getReadablePath(getParent(indexPrefix));

  # Argument 'pathnameD':
  pathnameD <- Arguments$getWritablePathname(pathnameD);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running BWA 'samse'");

  # Assert that input files are not overwritten
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameSAI));
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFQ));
##  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFA));

  res <- systemBWA("samse", "f"=shQuote(pathnameD), shQuote(indexPrefix), shQuote(pathnameSAI), shQuote(pathnameFQ), ..., verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # bwaSamse()


############################################################################
# HISTORY:
# 2014-03-10 [HB]
# o ROBUSTNESS: Now bwaSamse() uses shQuote() for all pathnames.
# 2012-09-24
# o Created.
############################################################################

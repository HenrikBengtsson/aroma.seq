###########################################################################/**
# @RdocDefault bwaSamse
#
# @title "Calls the BWA 'samse' command"
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
#   \item{pathnameFA}{The FASTA reference file.}
#   \item{pathnameD}{The destination pathname.}
#   \item{...}{Additional arguments specifying BWA 'samse' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   bwaIndex(pathnameFA)
#
#   pathnameSAI <- "bwaData/LambdaVirusExample/Generic/reads_1.sai";
#   pathnameFQ <- "fastqData/LambdaVirusExample/Generic/reads_1.fq";
#   pathnameD <- "bwaData/LambdaVirusExample/Generic/reads_1.sam";
#   bwaSamse(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ, pathnameFA=pathnameFA, pathnameD=pathnameD);
# }}
#
# @author
#*/###########################################################################
setMethodS3("bwaSamse", "default", function(pathnameSAI, pathnameFQ, pathnameFA, pathnameD, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameSAI':
  pathnameSAI <- Arguments$getReadablePathname(pathnameSAI);

  # Argument 'pathnameFQ':
  pathnameFQ <- Arguments$getReadablePathname(pathnameFQ);

  # Argument 'pathnameFA':
  pathnameFA <- Arguments$getReadablePathname(pathnameFA);

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
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFA));
 
  res <- systemBWA("samse", pathnameFA, pathnameSAI, pathnameFQ, "-f"=pathnameD, ..., verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # bwaSamse()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

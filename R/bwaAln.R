###########################################################################/**
# @RdocDefault bwaAln
#
# @title "BWA-backtrack alignment via 'bwa aln'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameFQ}{The FASTQ file to be aligned.}
#   \item{indexPrefix}{The pathname prefix to the BWA index files.}
#   \item{...}{Additional arguments specifying BWA 'aln' switches
#     passed to @see "systemBWA".}
#   \item{pathnameD}{The destination pathname.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
#   bwaIndex(pathnameFA)
#   indexPrefix <- bwaIndexPrefix(pathnameFA)
#   bwaAln("fastqData/LambdaVirusExample/Lambda_phage/reads_1.fq",
#          indexPrefix=indexPrefix,
#          pathnameD="fastqData/LambdaVirusExample/Lambda_phage/reads_1.sai")
# }}
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaAln", "default", function(pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


  verbose && enter(verbose, "Running BWA 'aln'");

  # Assert that input files are not overwritten
  .stop_if_not(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFQ));
##  .stop_if_not(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFA));

  res <- systemBWA("aln", "f"=shQuote(pathnameD), shQuote(indexPrefix), ..., shQuote(pathnameFQ), verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # bwaAln()

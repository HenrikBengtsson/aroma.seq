###########################################################################/**
# @RdocDefault bwaAln
#
# @title "Calls the BWA 'aln' command"
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
#   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   bwaIndex(pathnameFA, method="is")
#   indexPrefix <- bwaIndexPrefix(pathnameFA, method="is")
#   bwaAln("fastqData/LambdaVirusExample/Generic/reads_1.fq", indexPrefix=indexPrefix, pathnameD="fastqData/LambdaVirusExample/Generic/reads_1.sai")
# }}
#
# @author
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
  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFQ));
##  stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFA));
  
  res <- systemBWA("aln", "f"=pathnameD, indexPrefix, pathnameFQ, ..., verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # bwaAln()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

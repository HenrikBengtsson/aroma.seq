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
#   \item{filename, path}{The FASTQ file to be aligned.}
#   \item{pathnameD}{The destination pathname.}
#   \item{faPathname}{The FASTA reference file.}
#   \item{...}{Additional arguments specifying BWA 'aln' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{{\dontrun{
#   faPathname <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   bwaIndex(faPathname)
#   bwaAln("fastqData/LambdaVirusExample/Generic/reads_1.fq", pathnameD="fastqData/LambdaVirusExample/Generic/reads_1.sai", faPathname=faPathname)
# }}
#
# @author
#*/###########################################################################
setMethodS3("bwaAln", "default", function(filename, path=NULL, pathnameD, faPathname, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path);

  # Argument 'pathnameD':
  pathnameD <- Arguments$getWritablePathname(pathnameD);

  # Argument 'faPathname':
  faPathname <- Arguments$getReadablePathname(faPathname);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running BWA 'aln'");
  res <- systemBWA("aln", faPathname, pathname, stdout=pathnameD, ..., verbose=less(verbose, 10));
  verbose && exit(verbose);

  res;
}) # bwaAln()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

###########################################################################/**
# @RdocDefault bwaIndex
#
# @title "Calls the BWA index command"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{pathnameFA}{The FASTA file to be indexed.}
#   \item{...}{Additional arguments specifying BWA 'index' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   bwaIndex("annotationData/organisms/LambdaPhage/lambda_virus.fa")
# }}
#
# @author
#*/###########################################################################
setMethodS3("bwaIndex", "default", function(pathnameFA, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameFA':
  pathnameFA <- Arguments$getReadablePathname(pathnameFA);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running BWA index");
  res <- systemBWA("index", ..., pathnameFA, verbose=less(verbose, 10));
  verbose && exit(verbose);

  res;
}) # bwaIndex()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

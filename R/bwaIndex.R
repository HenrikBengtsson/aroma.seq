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
#   \item{filename, path}{The FASTA file to be indexed.}
#   \item{...}{Additional arguments specifying BWA 'index' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{{\dontrun{
#   bwaIndex("annotationData/organisms/LambdaPhage/lambda_virus.fa")
# }}
#
# @author
#*/###########################################################################
setMethodS3("bwaIndex", "default", function(filename, path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running BWA index");
  res <- systemBWA("index", ..., pathname, verbose=less(verbose, 10));
  verbose && exit(verbose);

  res;
}) # bwaIndex()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

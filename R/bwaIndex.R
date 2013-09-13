###########################################################################/**
# @RdocDefault bwaIndex
# @alias bwaIndex.FastaReferenceFile
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
#   \item{indexPrefix}{The prefix for the generated index files.}
#   \item{method}{Additional arguments passed to @see "bwaIndexPrefix".
#     Required if \code{indexPrefix == "*"}.}
#   \item{...}{Additional arguments specifying BWA 'index' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   bwaIndex("annotationData/organisms/LambdaPhage/lambda_virus.fa")
# }}
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaIndex", "default", function(pathnameFA, indexPrefix="*", method, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameFA':
  pathnameFA <- Arguments$getReadablePathname(pathnameFA);

  # Argument 'indexPrefix':
  if (identical(indexPrefix, "*")) {
    indexPrefix <- bwaIndexPrefix(pathnameFA, method=method);
  }
  if (!is.null(indexPrefix)) {
    path <- Arguments$getWritablePath(getParent(indexPrefix));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running BWA index");
  verbose && cat(verbose, "Index prefix: ", indexPrefix);
  res <- systemBWA("index", p=indexPrefix, ..., pathnameFA, verbose=less(verbose, 10));
  verbose && exit(verbose);

  res;
}) # bwaIndex()



############################################################################
# HISTORY:
# 2012-09-27
# o Added argument 'method' with will be passed to bwaIndexPrefix().
# 2012-09-24
# o Created.
############################################################################

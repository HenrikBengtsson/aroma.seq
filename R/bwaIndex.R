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
#   \item{...}{Additional arguments specifying BWA 'index' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   bwaIndex("annotationData/organisms/Lambda_phage/lambda_virus.fa")
# }}
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaIndex", "default", function(pathnameFA, indexPrefix="*", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameFA':
  pathnameFA <- Arguments$getReadablePathname(pathnameFA);

  # Argument 'indexPrefix':
  if (identical(indexPrefix, "*")) {
    indexPrefix <- bwaIndexPrefix(pathnameFA);
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
  res <- systemBWA("index", p=shQuote(indexPrefix), ..., shQuote(pathnameFA), verbose=less(verbose, 10));
  verbose && exit(verbose);

  res;
}) # bwaIndex()



############################################################################
# HISTORY:
# 2014-03-10 [HB]
# o ROBUSTNESS: Now bwaIndex() uses shQuote() for all pathnames.
# 2013-11-21
# o CLEANUP: bwaIndex() no longer has a 'method' argument as it is no
#   longer needed by bwaIndexPrefix().
# 2012-09-27
# o Added argument 'method' with will be passed to bwaIndexPrefix().
# 2012-09-24
# o Created.
############################################################################

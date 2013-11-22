###########################################################################/**
# @RdocDefault bwaIndexPrefix
#
# @title "Generates a prefix for the index files"
#
# \description{
#  @get "title" based on a FASTA pathname.
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameFA}{The FASTA file.}
#   \item{subdir}{The subdirectory relative to the FASTA file where to put
#     the BWA index files.}
#   \item{tags}{Tags added to the directory of the index set.}
#   \item{...}{Not used.}
# }
#
# \examples{
#   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   prefix <- bwaIndexPrefix(pathnameFA)
#   print(prefix)
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaIndexPrefix", "default", function(pathnameFA, subdir="bwa", tags="*", ...) {
  createIndexPrefix(pathnameFA, subdir=subdir, tags=tags, ...);
}) # bwaIndexPrefix()


############################################################################
# HISTORY:
# 2013-11-21
# o CLEANUP: bwaIndexPrefix() no longer adds a "methods" tags (i.e.
#   'bwtsw' or 'is'), because the generated index files are identical
#   regardless of method used.
# 2013-11-18
# o The default method for bwaIndexPrefix() is now 'bwtsw' (was 'is').
# 2012-09-25
# o Added argument 'method' and 'tags' to bwaIndexPrefix().
# 2012-09-24
# o Created.
############################################################################

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
#   \item{method}{The BWA algorithm used for generating the index set.}
#   \item{subdir}{The subdirectory relative to the FASTA file where to put
#     the BWA index files.}
#   \item{tags}{Tags added to the directory of the index set.}
#   \item{...}{Not used.}
# }
#
# \examples{
#   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   prefix <- bwaIndexPrefix(pathnameFA, method="is")
#   print(prefix)
# }
#
# @author
#*/###########################################################################
setMethodS3("bwaIndexPrefix", "default", function(pathnameFA, method=c("is", "bwtsw"), subdir="bwa", tags="*", ...) {
  createIndexPrefix(pathnameFA, subdir=subdir, tags=tags, asteriskTags=method, ...);
}) # bwaIndexPrefix()


############################################################################
# HISTORY:
# 2012-09-25
# o Added argument 'method' and 'tags' to bwaIndexPrefix().
# 2012-09-24
# o Created.
############################################################################

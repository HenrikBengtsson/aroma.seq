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
#   pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
#   prefix <- bwaIndexPrefix(pathnameFA)
#   print(prefix)
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaIndexPrefix", "default", function(pathnameFA, subdir="bwa", tags="*", ...) {
  createIndexPrefix(pathnameFA, subdir=subdir, tags=tags, ...)
}) # bwaIndexPrefix()

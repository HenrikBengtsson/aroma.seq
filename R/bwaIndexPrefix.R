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
#   \item{...}{Not used.}
# }
#
# \examples{
#   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
#   prefix <- bwaIndexPrefix(pathnameFA)
#   print(prefix)
# }
#
# @author
#*/###########################################################################
setMethodS3("bwaIndexPrefix", "default", function(pathnameFA, subdir="bwa", ...) {
  # Drop *.fa and *.fasta filename extensions
  prefix <- gsub("[.](fa|fasta)*$", "", pathnameFA, ignore.case=TRUE);
  path <- getParent(prefix);
  if (is.null(path)) {
    path <- subdir;
  } else {
    path <- file.path(path, subdir);
  }
  fullname <- basename(prefix);
  prefix <- file.path(path, fullname);
  prefix;
})


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################

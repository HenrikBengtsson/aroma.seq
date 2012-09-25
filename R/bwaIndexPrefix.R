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
  # Argument 'method':
  if (missing(method)) throw("Argument 'method' must be specified explicitly.");
  method <- match.arg(method);

  # Argument 'subdir':
  subdir <- Arguments$getCharacters(subdir);
  subdir <- do.call(file.path, as.list(subdir));

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);


  # Asterisk tags?
  keep <- (tags == "*");
  tags[keep] <- method;

  # Add tags to subdir
  subdir <- paste(c(subdir, tags), collapse=",");

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
}) # bwaIndexPrefix()


############################################################################
# HISTORY:
# 2012-09-25
# o Added argument 'method' and 'tags' to bwaIndexPrefix().
# 2012-09-24
# o Created.
############################################################################

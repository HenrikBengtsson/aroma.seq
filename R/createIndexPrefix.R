###########################################################################/**
# @RdocDefault createIndexPrefix
#
# @title "Generates a prefix for a set of index files"
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
#     the index files.}
#   \item{tags}{Tags added to the directory of the index set.}
#   \item{asteriskTags}{Tags to replace \code{"*"} in argument \code{tags}.}
#   \item{...}{Not used.}
# }
#
# \examples{
#   pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
#   prefix <- createIndexPrefix(pathnameFA)
#   print(prefix)
#   prefix <- createIndexPrefix(pathnameFA, tags="*,foo")
#   print(prefix)
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("createIndexPrefix", "default", function(pathnameFA, subdir=NULL, tags="*", asteriskTags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subdir':
  subdir <- Arguments$getCharacters(subdir);
  subdir <- do.call(file.path, as.list(subdir));

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Argument 'asteriskTags':
  if (!is.null(asteriskTags)) {
    asteriskTags <- Arguments$getTags(asteriskTags, collapse=",");
  }



  # Insert asterisk tags?
  if (length(tags) > 0L) {
    keep <- (tags == "*");
    if (any(keep)) {
      if (length(asteriskTags) > 0L) {
        tags[keep] <- asteriskTags;
      } else {
        tags[keep] <- NA;
      }
    }

    # Drop NA tags
    tags <- tags[!is.na(tags)];
  }

  # Add tags to subdir
  subdir <- paste(c(subdir, tags), collapse=",");
  subdir <- subdir[nchar(subdir) > 0L];

  # Drop *.fa and *.fasta filename extensions
  prefix <- gsub("[.](gz)*$", "", pathnameFA, ignore.case=TRUE);
  prefix <- gsub("[.](fa|fasta)*$", "", prefix, ignore.case=TRUE);
  path <- getParent(prefix);
  path <- Reduce(file.path, c(path, subdir));

  fullname <- basename(prefix);
  prefix <- file.path(path, fullname);

  prefix;
}) # createIndexPrefix()

###########################################################################/**
# @RdocClass BwaIndexSet
#
# @title "The BwaIndexSet class"
#
# \description{
#  @classhierarchy
#
#  An BwaIndexSet object represents a set of @see "BwaIndexFile":s.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "BwaIndexFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("BwaIndexSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "BwaIndexSet");
})


setMethodS3("as.character", "BwaIndexSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(this)));
  s <- c(s, sprintf("Complete: %s", isComplete(this)));

  class(s) <- class;
  s;
})


setMethodS3("getIndexPrefix", "BwaIndexSet", function(this, ...) {
  df <- getFile(this, 1L);
  getIndexPrefix(df, ...);
  path <- getPath(this);
  fullname <- getFullName(df);
  prefix <- file.path(path, fullname);
  prefix;
})


setMethodS3("isComplete", "BwaIndexSet", function(this, ...) {
  knownExts <- c("amb", "ann", "bwt", "pac", "sa");
  exts <- sapply(this, getExtension);
  missing <- setdiff(knownExts, exts);
  if (any(missing)) {
    return(FALSE);
  }
  TRUE;
})



############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################

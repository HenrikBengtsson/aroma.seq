###########################################################################/**
# @RdocClass AbstractIndexSet
#
# @title "The AbstractIndexSet class"
#
# \description{
#  @classhierarchy
#
#  An AbstractIndexSet object represents a set of @see "AbstractIndexFile":s.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AbstractIndexFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AbstractIndexSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "AbstractIndexSet");
})


setMethodS3("as.character", "AbstractIndexSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  class <- class(s);

  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(this)));
  s <- c(s, sprintf("Complete: %s", isComplete(this)));

  class(s) <- class;
  s;
}, protected=TRUE)


setMethodS3("byPrefix", "AbstractIndexSet", function(static, prefix, ...) {
  path <- getParent(prefix);
  byPath(static, path=path, ...);
}, static=TRUE)


setMethodS3("getIndexPrefix", "AbstractIndexSet", function(this, ...) {
  if (length(this) == 0L) return(as.character(NA));
  df <- getFile(this, 1L);
  getIndexPrefix(df, ...);
})


setMethodS3("isComplete", "AbstractIndexSet", abstract=TRUE);


############################################################################
# HISTORY:
# 2012-09-27
# o Created from BwaIndexSet.R
############################################################################

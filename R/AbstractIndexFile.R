###########################################################################/**
# @RdocClass AbstractIndexFile
#
# @title "The abstract AbstractIndexFile class"
#
# \description{
#  @classhierarchy
#
#  An AbstractIndexFile object represents a index file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \references{
#    ...
# }
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "AbstractIndexSet".
# }
#*/###########################################################################
setConstructorS3("AbstractIndexFile", function(...) {
  extend(GenericDataFile(...), "AbstractIndexFile");
})

setMethodS3("as.character", "AbstractIndexFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  class <- class(s);

  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(this)));

  class(s) <- class;
  s;
})


setMethodS3("getIndexPrefix", "AbstractIndexFile", function(this, ...) {
  path <- getPath(this);
  fullname <- getFullName(this);
  prefix <- file.path(path, fullname);
  prefix;
})


############################################################################
# HISTORY:
# 2012-09-27
# o Created from BwaIndexFile.R.
############################################################################

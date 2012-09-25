###########################################################################/**
# @RdocClass BwaIndexFile
#
# @title "The abstract BwaIndexFile class"
#
# \description{
#  @classhierarchy
#
#  A BwaIndexFile object represents a BWA index file.
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
#   @see "BwaIndexSet".
# }
#*/###########################################################################
setConstructorS3("BwaIndexFile", function(...) {
  extend(GenericDataFile(...), "BwaIndexFile");
})

setMethodS3("as.character", "BwaIndexFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(this)));

  class(s) <- class;
  s;
})


setMethodS3("getIndexPrefix", "BwaIndexFile", function(this, ...) {
  path <- getPath(this);
  fullname <- getFullName(this);
  prefix <- file.path(path, fullname);
  prefix;
})


############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################

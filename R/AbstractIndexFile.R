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
# @author "HB"
#
# \references{
#    ...
# }
#
# \seealso{
#   An object of this class is typically part of an
#   @see "AbstractIndexSet".
# }
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AbstractIndexFile", function(...) {
  extend(GenericDataFile(...), "AbstractIndexFile")
})

setMethodS3("as.character", "AbstractIndexFile", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(x)))
  s
}, protected=TRUE)


setMethodS3("getIndexPrefix", "AbstractIndexFile", function(this, ...) {
  path <- getPath(this)
  fullname <- getFullName(this)
  prefix <- file.path(path, fullname)
  prefix
})

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
#   \item{...}{Arguments passed to @see "AbstractIndexFile".}
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
#   @see "BwaIndexSet".
# }
#
# @keyword internal
#*/###########################################################################
setConstructorS3("BwaIndexFile", function(...) {
  extend(AbstractIndexFile(...), "BwaIndexFile")
})

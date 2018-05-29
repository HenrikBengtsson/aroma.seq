###########################################################################/**
# @RdocClass Bowtie2IndexFile
#
# @title "The abstract Bowtie2IndexFile class"
#
# \description{
#  @classhierarchy
#
#  A Bowtie2IndexFile object represents a Bowtie2 index file.
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
#   @see "Bowtie2IndexSet".
# }
#
# @keyword internal
#*/###########################################################################
setConstructorS3("Bowtie2IndexFile", function(...) {
  extend(AbstractIndexFile(...), "Bowtie2IndexFile")
})


setMethodS3("getIndexPrefix", "Bowtie2IndexFile", function(this, ...) {
  path <- getPath(this)
  fullname <- getFullName(this)
  fullname <- gsub("(|[.]rev)[.][0-9]$", "", fullname)
  prefix <- file.path(path, fullname)
  prefix
})

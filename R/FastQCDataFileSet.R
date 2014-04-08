###########################################################################/**
# @RdocClass FastQCDataFileSet
#
# @title "The FastQCDataFileSet class"
#
# \description{
#  @classhierarchy
#
#  An FastQCDataFileSet object represents a set of @see "FastQCDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "FastQCDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FastQCDataFileSet", function(files=NULL, ...) {
  extend(AromaSeqDataFileSet(files=files, ...), "FastQCDataFileSet")
})

setMethodS3("byPath", "FastQCDataFileSet", function(static, ..., recursive=TRUE, pattern="fastqc_data.txt$", verbose=FALSE) {
  res <- NextMethod("byPath", recursive=recursive, pattern=pattern);
  res;
}, static=TRUE)


############################################################################
# HISTORY:
# 2014-03-02
# o Created.
############################################################################

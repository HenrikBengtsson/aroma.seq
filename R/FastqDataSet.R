###########################################################################/**
# @RdocClass FastqDataSet
#
# @title "The FastqDataSet class"
#
# \description{
#  @classhierarchy
#
#  An FastqDataSet object represents a set of @see "FastqDataFile":s.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "FastqDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("FastqDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "FastqDataSet");
})


setMethodS3("validate", "FastqDataSet", function(this, ...) {
  NextMethod("validate", this, ...);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-06-28
# o Created.
############################################################################

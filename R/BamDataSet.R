###########################################################################/**
# @RdocClass BamDataSet
#
# @title "The BamDataSet class"
#
# \description{
#  @classhierarchy
#
#  An BamDataSet object represents a set of @see "BamDataFile":s.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "BamDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("BamDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "BamDataSet");
})


setMethodS3("validate", "BamDataSet", function(this, ...) {
  NextMethod("validate", this, ...);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-06-28
# o Created.
############################################################################

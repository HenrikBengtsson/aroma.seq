###########################################################################/**
# @RdocClass SamDataSet
#
# @title "The SamDataSet class"
#
# \description{
#  @classhierarchy
#
#  An SamDataSet object represents a set of @see "SamDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "SamDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SamDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "SamDataSet");
})


setMethodS3("validate", "SamDataSet", function(this, ...) {
  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getDepth", "SamDataSet", function(this, ...) {
  1L;
}, protected=TRUE);


setMethodS3("byPath", "SamDataSet", function(static, ..., pattern="[.](sam|SAM)$") {
  NextMethod("byPath", pattern=pattern);
}, static=TRUE)



############################################################################
# HISTORY:
# 2012-09-25
# o Created from BamDataSet.R.
############################################################################

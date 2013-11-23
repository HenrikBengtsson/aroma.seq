###########################################################################/**
# @RdocClass SamTransform
#
# @title "The SamTransform class"
#
# \description{
#  @classhierarchy
#
#  A SamTransform is an @see "AromaSeqTransform" that takes
#  @see "BamDataSet":s (or @see "SamDataSet":s) as input.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AromaSeqTransform".}
#  \item{.className}{A @character string specifying what class
#   of data sets to accept.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("SamTransform", function(..., .className="BamDataSet") {
  extend(AromaSeqTransform(..., .className=.className), "SamTransform")
}, abstract=TRUE)



setMethodS3("getRootPath", "SamTransform", function(this, ...) {
  # Use same root path as input data set, e.g. samData/ or bamData/
  ds <- getInputDataSet(this);
  path <- getPath(ds);
  path <- getParent(path, depth=2L);
  # Sanity check
  stopifnot(regexpr("Data$", path) != -1L);
  path;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-11-22
# o Created.
############################################################################

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
  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getDepth", "FastqDataSet", function(this, ...) {
  1L;
}, protected=TRUE);

setMethodS3("getDefaultSamReadGroup", "FastqDataSet", function(this, ...) {
  SamReadGroup();
})

setMethodS3("setSamReadGroup", "FastqDataSet", function(this, rg, ...) {
  # Argument 'rg':
  if (!is.null(rg)) {
    rg <- Arguments$getInstanceOf(rg, "SamReadGroup");
  }
  this$.rg <- rg;
  invisible(this);
})

setMethodS3("getSamReadGroup", "FastqDataSet", function(this, ...) {
  rg <- this$.rg;
  if (is.null(rg)) {
    rg <- getDefaultSamReadGroup(this, ...);
  }
  rg;
})


############################################################################
# HISTORY:
# 2012-06-28
# o Created.
############################################################################

##########################################################################/**
# @RdocClass BamIndexDataFile
#
# @title "The abstract BamIndexDataFile class"
#
# \description{
#  @classhierarchy
#
#  A BamIndexDataFile object represents a BAM index file.
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
# @keyword internal
#*/###########################################################################
setConstructorS3("BamIndexDataFile", function(...) {
  extend(GenericDataFile(...), "BamIndexDataFile");
})


setMethodS3("as.character", "BamIndexDataFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2012-10-02
# o Created.
############################################################################

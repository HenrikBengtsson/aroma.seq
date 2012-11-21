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
# @author
#*/###########################################################################
setConstructorS3("BamIndexDataFile", function(...) {
  extend(GenericDataFile(...), "BamIndexDataFile");
})


setMethodS3("as.character", "BamIndexDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  class <- class(s);

  class(s) <- class;
  s;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2012-10-02
# o Created.
############################################################################

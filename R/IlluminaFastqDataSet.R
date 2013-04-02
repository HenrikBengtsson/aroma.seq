###########################################################################/**
# @RdocClass IlluminaFastqDataSet
#
# @title "The IlluminaFastqDataSet class"
#
# \description{
#  @classhierarchy
#
#  An IlluminaFastqDataSet object represents a set of
#  @see "IlluminaFastqDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "FastqDataSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("IlluminaFastqDataSet", function(...) {
  extend(FastqDataSet(...), "IlluminaFastqDataSet");
})


############################################################################
# HISTORY:
# 2012-06-29
# o Created.
############################################################################

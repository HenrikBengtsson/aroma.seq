###########################################################################/**
# @RdocClass PairedEndFastqDataSet
#
# @title "The PairedEndFastqDataSet class"
#
# \description{
#  @classhierarchy
#
#  An PairedEndFastqDataSet object represents a paired-end @see "FastqDataSet".
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
setConstructorS3("PairedEndFastqDataSet", function(...) {
  extend(FastqDataSet(..., paired=TRUE), "PairedEndFastqDataSet")
})


setMethodS3("byPath", "PairedEndFastqDataSet", function(static, ...) {
  NextMethod("byPath", paired=TRUE);
}, protected=TRUE)

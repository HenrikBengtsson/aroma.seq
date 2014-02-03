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



############################################################################
# HISTORY:
# 2014-02-03
# o BUG FIX: byPath() for PairedEndFastqDataSet explicitly passed
#   arguments '...' to NextMethod(), which would cause them to be
#   duplicated in certain cases.
# 2013-11-04
# o A simple wrapper for paired-end FASTQ data sets.
# o Created.
############################################################################

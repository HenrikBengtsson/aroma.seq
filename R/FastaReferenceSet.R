###########################################################################/**
# @RdocClass FastaReferenceSet
#
# @title "The FastaReferenceSet class"
#
# \description{
#  @classhierarchy
#
#  An FastaReferenceSet object represents a set of @see "FastaReferenceFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "FastaReferenceFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "TT"
#*/###########################################################################
setConstructorS3("FastaReferenceSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "FastaReferenceSet");
})


setMethodS3("validate", "FastaReferenceSet", function(this, ...) {
  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getDepth", "FastaReferenceSet", function(this, ...) {
  1L;
}, protected=TRUE);


setMethodS3("byPath", "FastaReferenceSet", function(static, ..., pattern="[.](fa|fasta)$") {
  NextMethod("byPath", pattern=pattern);
}, static=TRUE)


############################################################################
# HISTORY:
# 2013-08-01 Created from SamDataSet.R
############################################################################

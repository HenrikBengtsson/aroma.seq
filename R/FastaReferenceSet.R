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

setMethodS3("getDepth", "FastaReferenceSet", function(this, ...) {
  1L;
}, protected=TRUE);


setMethodS3("byPath", "FastaReferenceSet", function(static, ..., pattern="[.](fa|fasta)(|[.]gz)$") {
  NextMethod("byPath", pattern=pattern);
}, static=TRUE)

setMethodS3("getOrganism", "FastaReferenceSet", function(this, ...) {
  aFile <- this[[1L]];
  getOrganism(aFile, ...);
})


############################################################################
# HISTORY:
# 2013-11-01
# o Now FastaReferenceSet$byPath() also finds gzip'ed FASTA files.
# 2013-08-01
# o Created from SamDataSet.R
############################################################################

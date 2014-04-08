###########################################################################/**
# @RdocClass SamDataFile
#
# @title "The abstract SamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A SamDataFile object represents a Sequence Alignment/Map (SAM) file [1].
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
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
# }
#
# \seealso{
#   An object of this class is typically part of an
#   @see "SamDataSet".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SamDataFile", function(...) {
  extend(AromaSeqDataFile(...), "SamDataFile");
})


setMethodS3("as.character", "SamDataFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s;
}, protected=TRUE)


setMethodS3("validate", "SamDataFile", validate.BamDataFile)


############################################################################
# HISTORY:
# 2013-11-16
# o Added validate() for SamDataFile.
# 2013-11-08
# o DOCUMENTATION: Added help on convertToBam() for SamDataFile.
# o Renamed to use convertToBam() for both SamDataFile and SamDataSet.
# 2012-09-25
# o Created from BamDataFile.R.
############################################################################

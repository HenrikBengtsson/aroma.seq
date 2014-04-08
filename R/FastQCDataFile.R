###########################################################################/**
# @RdocClass FastQCDataFile
#
# @title "The abstract FastQCDataFile class"
#
# \description{
#  @classhierarchy
#
#  A FastQCDataFile object represents a FastQC [1] data file.
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
# \references{
#  [1] Simon Andrews,
#      FastQC - A quality control tool for high throughput sequence data,
#      March 2014.
#      \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
# }
#
# \seealso{
#   An object of this class is typically part of a @see "FastQCDataFileSet".
# }
#*/###########################################################################
setConstructorS3("FastQCDataFile", function(...) {
  extend(AromaSeqDataFile(...), "FastQCDataFile");
})

setMethodS3("as.character", "FastQCDataFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s;
}, protected=TRUE)

setMethodS3("getSampleName", "FastQCDataFile", function(this, ...) {
  name <- getPath(this);
  name <- basename(name);
  name <- gsub("_fastqc$", "", name);
  name;
}, protected=TRUE)

setMethodS3("getDefaultFullName", "FastQCDataFile", function(this, ...) {
  getSampleName(this, ...);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2014-03-02
# o Created.
############################################################################

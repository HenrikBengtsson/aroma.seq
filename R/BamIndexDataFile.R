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

  s <- NextMethod("as.character", ...);
  class <- class(s);

  class(s) <- class;
  s;
})


setMethodS3("extractStats", "BamIndexDataFile", function(this, ..., force=FALSE) {
  stats <- this$.stats;

  if (force || is.null(stats)) {
    pathname <- getPathname(this);
    bfr <- systemSamtools("idxstats", pathname, stdout=TRUE, ...);
    stats <- strsplit(bfr, split="\t", fixed=TRUE);
    seqName <- sapply(stats, FUN=.subset, 1L);
    seqLength <- sapply(stats, FUN=.subset, 2L);
    seqLength <- Arguments$getIntegers(seqLength);
    countMapped <- sapply(stats, FUN=.subset, 3L);
    countMapped <- Arguments$getIntegers(countMapped);
    countUnmapped <- sapply(stats, FUN=.subset, 4L);
    countUnmapped <- Arguments$getIntegers(countUnmapped);
    data <- data.frame(length=seqLength, mapped=countMapped, unmapped=countUnmapped);
    rownames(data) <- seqName;
    stats <- data;
    this$.stats <- stats;
  }

  data;
})


############################################################################
# HISTORY:
# 2012-10-02
# o Created.
############################################################################

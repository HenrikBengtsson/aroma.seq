
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
  this <- extend(FastqDataSet(...), "IlluminaFastqDataSet")
  validate(this)
})


setMethodS3("as.character", "IlluminaFastqDataSet", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, sprintf("Platform: %s", htable(sapply(x, FUN=getPlatform))))
  s <- c(s, sprintf("File format version: %s", htable(sapply(x, FUN=getFileVersion))))
  s
}, protected=TRUE)


setMethodS3("validate", "IlluminaFastqDataSet", function(this, ...) {
  NextMethod("validate")
  for (ii in seq_along(this)) {
    fq <- this[[ii]]
    if (!inherits(fq, "IlluminaFastqDataFile")) {
      throw(sprintf("Invalid %s: Detected non-IlluminaFastqDataFile (of class %s): %s", class(this)[1], class(fq)[1], getPathname(fq)))
    }
    if (is.na(getFileVersion(fq))) {
      throw(sprintf("Invalid %s: Detected %s with unknown file format version: %s", class(this)[1], class(fq)[1], getPathname(fq)))
    }
  }
  invisible(this)
}, protected=TRUE)


############################################################################
# HISTORY:
# 2012-06-29
# o Created.
############################################################################

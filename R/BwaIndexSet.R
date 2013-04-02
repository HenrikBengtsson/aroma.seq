###########################################################################/**
# @RdocClass BwaIndexSet
#
# @title "The BwaIndexSet class"
#
# \description{
#  @classhierarchy
#
#  An BwaIndexSet object represents a set of @see "BwaIndexFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AbstractIndexSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("BwaIndexSet", function(...) {
  extend(AbstractIndexSet(...), "BwaIndexSet");
})

setMethodS3("isComplete", "BwaIndexSet", function(this, ...) {
  knownExts <- c("amb", "ann", "bwt", "pac", "sa");
  if (length(this) < length(knownExts)) return(FALSE);

  exts <- sapply(this, getExtension);
  missing <- setdiff(knownExts, exts);
  if (any(missing)) {
    return(FALSE);
  }

  TRUE;
})



############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################

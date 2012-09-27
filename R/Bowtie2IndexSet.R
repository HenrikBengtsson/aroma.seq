###########################################################################/**
# @RdocClass Bowtie2IndexSet
#
# @title "The Bowtie2IndexSet class"
#
# \description{
#  @classhierarchy
#
#  An Bowtie2IndexSet object represents a set of @see "Bowtie2IndexFile":s.
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
# @author
#*/###########################################################################
setConstructorS3("Bowtie2IndexSet", function(...) {
  extend(AbstractIndexSet(...), "Bowtie2IndexSet");
})

setMethodS3("isComplete", "Bowtie2IndexSet", function(this, ...) {
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

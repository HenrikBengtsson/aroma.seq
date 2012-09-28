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
  knownExts <- c("1", "2", "3", "4", "rev.1", "rev.2");
  knownExts <- sprintf("%s.bt2", knownExts);
  if (length(this) < length(knownExts)) return(FALSE);

  filenames <- sapply(this, getFilename);
  patterns <- sprintf("[.]%s$", knownExts);
  idxs <- sapply(patterns, FUN=grep, filenames);
  ns <- sapply(idxs, FUN=length);

  if(any(ns == 0L)) return(FALSE);

  TRUE;
})



############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################

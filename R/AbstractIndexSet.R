###########################################################################/**
# @RdocClass AbstractIndexSet
#
# @title "The AbstractIndexSet class"
#
# \description{
#  @classhierarchy
#
#  An AbstractIndexSet object represents a set of @see "AbstractIndexFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AbstractIndexFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
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
setConstructorS3("AbstractIndexSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "AbstractIndexSet");
})


setMethodS3("as.character", "AbstractIndexSet", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s <- c(s, sprintf("Index prefix: %s", getIndexPrefix(this)));
  s <- c(s, sprintf("Organism: %s", getOrganism(this)));
  s <- c(s, sprintf("Complete: %s", isComplete(this)));
  s;
}, protected=TRUE)


setMethodS3("getOrganism", "AbstractIndexSet", function(this, ...) {
  path <- getPath(this);
  path <- dirname(path);
  organism <- basename(path);
  organism;
})


setMethodS3("byPrefix", "AbstractIndexSet", function(static, prefix, ...) {
  path <- getParent(prefix);
  pattern <- sprintf("%s", basename(prefix));
  byPath(static, path=path, pattern=pattern, ...);
}, static=TRUE)


setMethodS3("getIndexPrefix", "AbstractIndexSet", function(this, ...) {
  df <- getOneFile(this);
  getIndexPrefix(df, ...);
})


setMethodS3("getFastaReferenceFile", "AbstractIndexSet", function(this, ...) {
  organism <- getOrganism(this);
  prefix <- getIndexPrefix(this);
  fullname <- basename(prefix);
  fa <- FastaReferenceFile$byOrganism(organism, prefix=fullname);

  # Assert compatibility
  isCompatibleWith(fa, this, mustWork=TRUE);

  fa;
})


setMethodS3("isComplete", "AbstractIndexSet", abstract=TRUE);


############################################################################
# HISTORY:
# 2014-08-23
# o ROBUSTNESS: Now getFastaReferenceFile() now asserts compatibility.
# 2014-01-18
# o Added getFastaReferenceFile() for AbstractIndexSet.
# 2013-11-17
# o BUG FIX: BwaIndexSet$byPrefix(prefix) would find any BWA index set
#   in directory dirname(prefix) without matching filenames of the set
#   to basename(prefix).
# 2013-11-10
# o Added getOrganism().
# 2012-09-27
# o Created from BwaIndexSet.R
############################################################################

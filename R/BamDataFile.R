###########################################################################/**
# @RdocClass BamDataFile
#
# @title "The abstract BamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A BamDataFile object represents a BAM file.
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
#
# \references{
#    ...
# }
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "BamDataSet".
# }
#*/###########################################################################
setConstructorS3("BamDataFile", function(...) {
  extend(GenericDataFile(...), "BamDataFile");
})


setMethodS3("as.character", "BamDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  n <- nbrOfTargets(this);
  s <- c(s, sprintf("Number of targets: %s", n));
  len <- getTotalTargetLength(this);
  s <- c(s, sprintf("Total target length: %.3gMb (%.0f bases)", len/1e9, len));
  names <- getTargetNames(this);
  s <- c(s, sprintf("Targets: [%d] %s", n, hpaste(names)));
  s <- c(s, sprintf("Has index file (*.bai): %s", hasIndex(this)));
  s <- c(s, sprintf("Is sorted: %s", isSorted(this)));

  class(s) <- class;
  s;
})

setMethodS3("hasIndex", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  pathnameI <- sprintf("%s.bai", pathname);
  isFile(pathnameI);
})

# \details{
#   BAM headers typically contain an \code{"@HD VN:1.0 SO:<value>"} entry,
#   where \code{<value>} indicates whether the aligned reads are sorted 
#   or not.  Unfortunately, this entry is neither enforced nor has it to
#   be correct [1,2].
#
#   Instead, we consider a BAM file to be sorted if and only if it has
#   an index file.  The rationale is that it is not possible to index
#   a BAM file unless it is sorted first.
# }
#
# \references{
#   [1] Question: is my BAM file sorted?, Biostar, 2011,
#   \url{http://www.biostars.org/post/show/5256/is-my-bam-file-sorted/}\cr
#   [2] Asking for suggestiona on samtools bug fixing, SEQanswers, 2010,
#   \url{http://seqanswers.com/forums/showthread.php?t=3739}\cr
# }
setMethodS3("isSorted", "BamDataFile", function(this, ...) {
  isTRUE(hasIndex(this));
})


setMethodS3("sort", "BamDataFile", function(this, force=FALSE, ...) {
  # Nothing todo?
  if (!force && isSorted(this)) {
    return(this);
  }

  throw("Not yet implemented!");
})


setMethodS3("nbrOfSeqs", "BamDataFile", function(this, ...) {
  nbrOfTargets(this);
})

setMethodS3("getTargets", "BamDataFile", function(this, ...) {
  hdr <- getHeader(this);
  targets <- hdr$targets;
  targets;  
})

setMethodS3("nbrOfTargets", "BamDataFile", function(this, ...) {
  length(getTargets(this));
})

setMethodS3("getTargetNames", "BamDataFile", function(this, ...) {
  names(getTargets(this));
})

setMethodS3("getTargetLengths", "BamDataFile", function(this, ...) {
  getTargets(this);
})

setMethodS3("getTotalTargetLength", "BamDataFile", function(this, ...) {
  sum(as.numeric(getTargets(this)));
})

setMethodS3("getHeader", "BamDataFile", function(this, force=FALSE, ...) {
  header <- this$.header;
  if (force || is.null(header)) {
    header <- readHeader(this, ...);
    this$.header <- header;
  }
  header;
})

setMethodS3("readHeader", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  bf <- Rsamtools::BamFile(pathname);
  hdr <- Rsamtools::scanBamHeader(bf);
  hdr;
}, private=TRUE)


############################################################################
# HISTORY:
# 2012-06-28
# o Created.
############################################################################

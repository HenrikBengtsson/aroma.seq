###########################################################################/**
# @RdocClass SamDataFile
#
# @title "The abstract SamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A SamDataFile object represents a SAM file.
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
#   @see "SamDataSet".
# }
#*/###########################################################################
setConstructorS3("SamDataFile", function(...) {
  extend(GenericDataFile(...), "SamDataFile");
})


setMethodS3("as.character", "SamDataFile", function(x, ...) {
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

setMethodS3("hasIndex", "SamDataFile", function(this, ...) {
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
setMethodS3("isSorted", "SamDataFile", function(this, ...) {
  isTRUE(hasIndex(this));
})


# Argument '...' must be 2nd to match the generic base::sort() function.
setMethodS3("sort", "SamDataFile", function(x, ..., force=FALSE) {
  # To please R CMD check
  this <- x;

  # Nothing todo?
  if (!force && isSorted(this)) {
    return(this);
  }

  throw("Not yet implemented!");
})


setMethodS3("nbrOfSeqs", "SamDataFile", function(this, ...) {
  nbrOfTargets(this);
})

setMethodS3("getTargets", "SamDataFile", function(this, ...) {
  hdr <- getHeader(this);
  targets <- hdr$targets;
  targets;  
})

setMethodS3("nbrOfTargets", "SamDataFile", function(this, ...) {
  length(getTargets(this));
})

setMethodS3("getTargetNames", "SamDataFile", function(this, ...) {
  names(getTargets(this));
})

setMethodS3("getTargetLengths", "SamDataFile", function(this, ...) {
  getTargets(this);
})

setMethodS3("getTotalTargetLength", "SamDataFile", function(this, ...) {
  sum(as.numeric(getTargets(this)));
})

setMethodS3("getHeader", "SamDataFile", function(this, force=FALSE, ...) {
  header <- this$.header;
  if (force || is.null(header)) {
    header <- readHeader(this, ...);
    this$.header <- header;
  }
  header;
})

setMethodS3("readHeader", "SamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  bf <- Rsamtools::SamFile(pathname);
  hdr <- Rsamtools::scanSamHeader(bf);
  hdr;
}, private=TRUE)


############################################################################
# HISTORY:
# 2012-09-25
# o Created from BamDataFile.R.
############################################################################

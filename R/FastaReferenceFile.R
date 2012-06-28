###########################################################################/**
# @RdocClass FastaReferenceFile
#
# @title "The FastaReferenceFile class"
#
# \description{
#  @classhierarchy
#
#  A FastaReferenceFile object represents a FASTA reference file.
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
# \seealso{
#   ...
# }
#*/###########################################################################
setConstructorS3("FastaReferenceFile", function(...) {
  extend(GenericDataFile(...), "FastaReferenceFile",
    .seqLengths=NULL
  );
})

setMethodS3("as.character", "FastaReferenceFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  n <- nbrOfSeqs(this);
  s <- c(s, sprintf("Total sequence length: %d", getTotalSeqLengths(this)));
  s <- c(s, sprintf("Number of sequences: %s", n));
  s <- c(s, sprintf("Sequence names: [%d] %s", n, hpaste(getSeqNames(this))));

  class(s) <- class;
  s;
})


setMethodS3("getSeqLengths", "FastaReferenceFile", function(this, force=FALSE, ...) {
  seqLengths <- this$.seqLengths;
  if (force || is.null(seqLengths)) {
    seqLengths <- readSeqLengths(this, ...);
    this$.seqLengths <- seqLengths;
  }
  seqLengths;
})

setMethodS3("getTotalSeqLengths", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  sum(seqLengths);
})

setMethodS3("getSeqNames", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  names(seqLengths);
})

setMethodS3("nbrOfSeqs", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  length(seqLengths);
})


# \seealso{
#   Internally, \code{fasta.info()} of \pkg{Biostrings} is used.
# }
setMethodS3("readSeqLengths", "FastaReferenceFile", function(this, ...) {
  pathname <- getPathname(this);
  seqLengths <- Biostrings::fasta.info(pathname);
  seqLengths;
}, private=TRUE)


############################################################################
# HISTORY:
# 2012-06-28
# o Created.
############################################################################

###########################################################################/**
# @RdocClass FastaReferenceIndexFile
#
# @title "The FastaReferenceIndexFile class"
#
# \description{
#  @classhierarchy
#
#  A FastaReferenceIndexFile object represents a
#  FASTA reference index file (FAI).
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
# \seealso{
#   [1] \url{http://www.wikipedia.org/wiki/FASTA_format}
# }
#*/###########################################################################
setConstructorS3("FastaReferenceIndexFile", function(...) {
  extend(GenericDataFile(...), "FastaReferenceIndexFile")
})

setMethodS3("as.character", "FastaReferenceIndexFile", function(x, ...) {
  this <- x
  s <- NextMethod("as.character")
  n <- nbrOfSeqs(this)
  s <- c(s, sprintf("Total sequence length: %s", pi3(getTotalSeqLengths(this))))
  s <- c(s, sprintf("Number of sequences: %d", n))
  s <- c(s, sprintf("Sequence names: [%d] %s", n, hpaste(getSeqNames(this))))
  s
}, protected=TRUE)


setMethodS3("getDefaultFullName", "FastaReferenceIndexFile", function(this, ...) {
  name <- NextMethod("getDefaultFullName")
  name <- gsub("[.](fa|fasta)[.](fai)$", "", name, ignore.case=TRUE)
  name
}, protected=TRUE)


setMethodS3("getSeqLengths", "FastaReferenceIndexFile", function(this, force=FALSE, ...) {
  readSeqLengths(this, ...)
})

setMethodS3("getTotalSeqLengths", "FastaReferenceIndexFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...)
  if (is.null(seqLengths)) return(NA_integer_)
  res <- sum(as.numeric(seqLengths))
  if (res < .Machine$integer.max) res <- as.integer(res)
  res
})

setMethodS3("getSeqNames", "FastaReferenceIndexFile", function(this, ...) {
  names(getSeqLengths(this, ...))
})

setMethodS3("nbrOfSeqs", "FastaReferenceIndexFile", function(this, ...) {
  length(getSeqLengths(this, ...))
})


setMethodS3("readDataFrame", "FastaReferenceIndexFile", function(this, ...) {
  pathname <- getPathname(this)
  data <- read.table(pathname, sep="\t", colClasses=c("character", rep("integer", 4L)))
  colnames(data) <- c("sequence", "length", "fileOffset", "lengthPerEntry", "bytesPerEntry")
  data
}, private=TRUE)

setMethodS3("readSeqLengths", "FastaReferenceIndexFile", function(this, force=FALSE, ...) {
  data <- readDataFrame(this, ...)
  lens <- data$length
  names(lens) <- data$sequence
  lens
}, private=TRUE)

############################################################################
# HISTORY:
# 2015-05-01
# o Created.
############################################################################

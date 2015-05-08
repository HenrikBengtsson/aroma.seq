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
  extend(GenericDataFile(...), c("FastaReferenceIndexFile", uses("SequenceContigsInterface")))
})

setMethodS3("as.character", "FastaReferenceIndexFile", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, getSeqGenericSummary(x, ...))
  s
}, protected=TRUE)


setMethodS3("getDefaultFullName", "FastaReferenceIndexFile", function(this, ...) {
  name <- NextMethod("getDefaultFullName")
  name <- gsub("[.](fa|fasta)[.](fai)$", "", name, ignore.case=TRUE)
  name
}, protected=TRUE)


setMethodS3("readSeqLengths", "FastaReferenceIndexFile", function(this, force=FALSE, ...) {
  data <- readDataFrame(this, ...)
  lens <- data$size
  names(lens) <- data$contig
  lens
}, private=TRUE)

setMethodS3("getSeqLengths", "FastaReferenceIndexFile", function(this, force=FALSE, ...) {
  readSeqLengths(this, ...)
})

setMethodS3("readDataFrame", "FastaReferenceIndexFile", function(this, ...) {
  pathname <- getPathname(this)
  data <- read.table(pathname, sep="\t", colClasses=c("character", rep("integer", 4L)))
  ## Column names as in GATK (suggested in thread 'How can I prepare a
  ## FASTA file to use as reference?' last edited Sept 2013
  ## [http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference]
  colnames(data) <- c("contig", "size", "location", "basesPerLine", "bytesPerLine")
  data
}, private=TRUE)


############################################################################
# HISTORY:
# 2015-05-01
# o Created.
############################################################################

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
# @author "HB"
#
# \references{
#  [1] Trapnell et al. \emph{Differential gene and transcript expression
#      analysis of RNA-seq experiments with TopHat and Cufflinks}.
#      Nat Protoc, 2012.\cr
# }
#
# @keyword internal
#*/###########################################################################
setConstructorS3("Bowtie2IndexSet", function(...) {
  extend(AbstractIndexSet(...), "Bowtie2IndexSet");
})


setMethodS3("isComplete", "Bowtie2IndexSet", function(this, ...) {
  # TODO: The set of index files generated may vary with
  # bowtie2-index options. /HB 2012-09-27
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


setMethodS3("getSummary", "Bowtie2IndexSet", function(this, ...) {
  stopifnot(isComplete(this));
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  prefix <- getIndexPrefix(this);
  tempfile
  bfr <- system2("bowtie2-inspect", args=c("--summary", prefix), stdout=TRUE);
  keys <- gsub("\\t.*", "", bfr);
  values <- strsplit(bfr, split="\t", fixed=TRUE);
  values <- lapply(values, FUN="[", -1L);
  names(values) <- keys;
  values;
})


setMethodS3("getSeqNames", "Bowtie2IndexSet", function(this, ...) {
  stopifnot(isComplete(this));
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  prefix <- getIndexPrefix(this);
  names <- system2("bowtie2-inspect", args=c("--names", prefix), stdout=TRUE);

  names <- trim(names);

  # From http://en.wikipedia.org/wiki/FASTA_format:
  # "The **word** following the ">" symbol is the identifier of the
  # sequence, and the rest of the line is the description [...]", e.g.
  # >I dna:chromosome chromosome:EF4:I:1:230218:1 REF
  # => ID/name: 'I'
  # => Description: 'dna:chromosome chromosome:EF4:I:1:230218:1 REF'
  names <- gsub(" .*", "", names);

  names;
})


setMethodS3("getSequenceNames", "Bowtie2IndexSet", function(this, ...) {
  .Deprecated("getSeqNames")
  getSeqNames(this, ...)
}, protected=TRUE, deprecated=TRUE)


############################################################################
# HISTORY:
# 2014-08-23
# o ROBUSTNESS: Now buildBowtie2IndexSet() asserts that the returned
#   index set is compatible with the FASTA file.
# 2014-08-11
# o BUG FIX: getSeqNames() for Bowtie2IndexSet would return the
#   sequence description in addition to the ID as part of the name.
# 2014-07-24
# o CONSISTENCY: Renamed getSequenceNames() to getSeqNames() for
#   Bowtie2IndexSet.  Deprecated the old version.
# 2014-04-10
# o BUG FIX: getSequenceNames() for Bowtie2IndexSet returned nothing.
# 2012-09-27
# o Added getSummary() and getSequenceNames() for Bowtie2IndexSet, which
#   utilizes 'bowtie2-inspect' executable.
# 2012-09-27
# o Created from BwaIndexSet.R.
############################################################################

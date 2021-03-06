###########################################################################/**
# @RdocClass Bowtie2IndexSet
# @aliasmethod getDefaultFilePatterns
# @aliasmethod readSeqLengths
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
  extend(AbstractIndexSet(...), c("Bowtie2IndexSet", uses("SequenceContigsInterface")))
})


setMethodS3("as.character", "Bowtie2IndexSet", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, getSeqGenericSummary(x, ...))
  s
})


setMethodS3("getDefaultFilePatterns", "Bowtie2IndexSet", function(static, prefix, ...) {
  knownExts <- c("1", "2", "3", "4", "rev.1", "rev.2")
  sprintf("%s[.](%s)[.]bt2$", basename(prefix), paste(knownExts, collapse="|"))
}, static=TRUE, protected=TRUE)


setMethodS3("isComplete", "Bowtie2IndexSet", function(this, ...) {
  # TODO: The set of index files generated may vary with
  # bowtie2-index options. /HB 2012-09-27
  knownExts <- c("1", "2", "3", "4", "rev.1", "rev.2")
  knownExts <- sprintf("%s.bt2", knownExts)
  if (length(this) < length(knownExts)) return(FALSE)

  filenames <- sapply(this, getFilename)
  patterns <- sprintf("[.]%s$", knownExts)
  idxs <- sapply(patterns, FUN=grep, filenames)
  ns <- sapply(idxs, FUN=length)

  if(any(ns == 0L)) return(FALSE)

  TRUE
})


setMethodS3("getSummary", "Bowtie2IndexSet", function(this, ...) {
  .stop_if_not(isComplete(this))
  .stop_if_not(isCapableOf(aroma.seq, "bowtie2"))
  prefix <- getIndexPrefix(this)
  tempfile
  bfr <- system2("bowtie2-inspect", args=c("--summary", prefix), stdout=TRUE)
  keys <- gsub("\\t.*", "", bfr)
  values <- strsplit(bfr, split="\t", fixed=TRUE)
  values <- lapply(values, FUN="[", -1L)
  names(values) <- keys
  values
})


setMethodS3("readSeqLengths", "Bowtie2IndexSet", function(this, force=FALSE, ...) {
  .stop_if_not(isComplete(this))
  .stop_if_not(isCapableOf(aroma.seq, "bowtie2"))

  prefix <- getIndexPrefix(this)

  # Check for cached results
  dirs <- c("aroma.seq", getOrganism(this))
  key <- list(method="readSeqLengths", class=class(this), prefix=prefix)
  lens <- loadCache(key=key, dirs=dirs)
  if (!force && !is.null(lens)) {
    return(lens)
  }

  names <- system2("bowtie2-inspect", args=c("--names", prefix), stdout=TRUE)
  names <- trim(names)

  lens <- rep(NA_integer_, times=length(names))
  names(lens) <- names

  # Cache
  saveCache(lens, key=key, dirs=dirs)

  lens
}, protected=TRUE)


setMethodS3("getSeqLengths", "Bowtie2IndexSet", function(this, force=FALSE, ...) {
  lens <- this$.seqLengths
  if (force || is.null(lens)) {
    lens <- readSeqLengths(this, ...)
    this$.seqLengths <- lens
  }
  lens
})


setMethodS3("isTopHat2IndexSet", "Bowtie2IndexSet", function(this, ...) {
  # AD HOC
  path <- getPath(this)
  grepl("tophat2", path, fixed=TRUE)
}, protected=TRUE)

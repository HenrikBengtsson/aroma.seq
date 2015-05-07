###########################################################################/**
# @RdocClass SequenceContigsInterface
#
# @title "The SequenceContigsInterface class"
#
# \description{
#  @classhierarchy
#
#  An @see "R.oo::Interface" for methods for with sequence contigs.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SequenceContigsInterface", function(...) {
  extend(Interface(...), "SequenceContigsInterface")
})


setMethodS3("getSeqGenericSummary", "SequenceContigsInterface", function(x, ...) {
  n <- nbrOfSeqs(x)
  s <- sprintf("Number of sequence contigs: %d", n)
  if (n > 0L) {
    names <- getSeqNames(x, unique=TRUE)
    dups <- hasDuplicatedSeqs(x)
    if (dups) {
      nu <- nbrOfSeqs(x, unique=TRUE)
      msg <- sprintf("Unique sequence names: [%d] %s", nu, hpaste(names))
    } else {
      msg <- sprintf("Sequence names: [%d] %s", n, hpaste(names))
    }
    s <- c(s, msg)
    if (hasSeqLengths(x)) {
      lens <- getSeqLengths(x)
      total <- getTotalSeqLength(x)
      s <- c(s, sprintf("Sequence lengths (bp): [%d] %s", n, hpaste(pi3(lens), collapse="; ")))
      s <- c(s, sprintf("Total sequence length (bp): %s", pi3(total)))
    }
    if (n > 1L) {
      scores <- getSeqOrdering(x, rank=TRUE, as="humanreadable")
      scores <- paste(scores, collapse=", ")
      s <- c(s, sprintf("Ordering of sequence names (scores): %s", scores))
    }
  }
  GenericSummary(s)
}, protected=TRUE)

setMethodS3("getSeqOrdering", "SequenceContigsInterface", function(this, ...) {
  lens <- getSeqLengths(this)
  typeOfSequenceOrdering(lens, ...)
}, protected=TRUE)

setMethodS3("getSeqLengths", "SequenceContigsInterface", abstract=TRUE)

setMethodS3("hasSeqLengths", "SequenceContigsInterface", function(this, ...) {
  lens <- getSeqLengths(this, ...)
  !all(is.na(lens))
})

setMethodS3("getTotalSeqLength", "SequenceContigsInterface", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...)
  if (is.null(seqLengths)) return(NA_integer_)
  res <- sum(as.numeric(seqLengths), na.rm=TRUE)
  if (res < .Machine$integer.max) res <- as.integer(res)
  res
})


setMethodS3("cleanSeqNames", "SequenceContigsInterface", function(this, names, ...) {
  # From http://en.wikipedia.org/wiki/FASTA_format:
  # "The **word** following the ">" symbol is the identifier of the
  # sequence, and the rest of the line is the description [...]", e.g.
  # >I dna:chromosome chromosome:EF4:I:1:230218:1 REF
  # => ID/name: 'I'
  # => Description: 'dna:chromosome chromosome:EF4:I:1:230218:1 REF'
  names <- gsub(" .*", "", names)
  names
}, protected=TRUE)


setMethodS3("getSeqNames", "SequenceContigsInterface", function(this, clean=TRUE, unique=FALSE, ...) {
  seqLengths <- getSeqLengths(this, ...)
  if (is.null(seqLengths)) return(NA_character_)
  names <- names(seqLengths)
  if (clean) names <- cleanSeqNames(this, names)
  if (unique) names <- unique(names)
  names
})

setMethodS3("nbrOfSeqs", "SequenceContigsInterface", function(this, ...) {
  names <- getSeqNames(this, ...)
  if (length(names) == 1L && is.na(names)) return(NA_integer_)
  length(names)
})

setMethodS3("hasDuplicatedSeqs", "SequenceContigsInterface", function(this, ...) {
  anyDuplicated(getSeqNames(this, unique=FALSE, ...))
}, protected=TRUE)

setMethodS3("isCompatibleWithBySeqLengths", "SequenceContigsInterface", function(this, other, ...) {
  if (!inherits(other, "SequenceContigsInterface")) {
    res <- FALSE
    attr(res, "reason") <- sprintf("The 'other' object is not an SequenceContigsInterface class: %s", class(other)[1L])
    return(res)
  }

  lens <- getSeqLengths(this)
  lensO <- getSeqLengths(other)
  idxs <- match(lensO, lens)
  res <- all(diff(idxs) > 0, na.rm=TRUE)
  if (!res) {
    attr(res, "reason") <- "The ordering of sequence lengths does not match."
  }

  res
}, protected=TRUE)


setMethodS3("isCompatibleWithBySeqNames", "SequenceContigsInterface", function(this, other, unique=TRUE, ...) {
  if (!inherits(other, "SequenceContigsInterface")) {
    res <- FALSE
    attr(res, "reason") <- sprintf("The 'other' object is not an SequenceContigsInterface class: %s", class(other)[1L])
    return(res)
  }

  names <- getSeqNames(this, unique=unique)
  namesO <- getSeqNames(other, unique=unique)
  idxs <- match(namesO, names)
  res <- all(diff(idxs) > 0, na.rm=TRUE)
  if (!res) {
    attr(res, "reason") <- "The ordering of sequence names does not match."
  }

  res
}, protected=TRUE)


setMethodS3("isCompatibleWithBySeqs", "SequenceContigsInterface", function(this, other, ...) {
  res <- isCompatibleWithBySeqNames(this, other, ...)
  if (!res) return(res)

  if (hasSeqLengths(this)) {
    res <- isCompatibleWithBySeqLengths(this, other, ...)
    if (!res) return(res)
  }

  res
}, protected=TRUE)



############################################################################
# HISTORY:
# 2015-05-06
# o Extracted from FastaReferenceFile.R.
############################################################################

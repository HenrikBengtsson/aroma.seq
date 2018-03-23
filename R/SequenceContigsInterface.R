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
  if (!is.na(n) && n > 0L) {
    dups <- hasDuplicatedSeqs(x)
    names <- getSeqNames(x, unique=dups)

    if (dups) {
      nu <- length(names)
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
      scores <- getSeqOrdering(x, unique=dups, rank=TRUE, as="humanreadable")
      scores <- paste(scores, collapse=", ")
      if (dups) {
        msg <- sprintf("Ordering of unique sequence contigs (scores): %s", scores)
      } else {
        msg <- sprintf("Ordering of sequence contigs (scores): %s", scores)
      }
      s <- c(s, msg)
    }
  }
  GenericSummary(s)
}, protected=TRUE)

setMethodS3("getSeqOrdering", "SequenceContigsInterface", function(this, ...) {
  lens <- getSeqLengths(this)
  if (!all(is.na(lens))) {
    seqOrder <- typeOfSequenceOrdering(lens, ...)
  } else {
    names <- names(lens)
    names <- cleanSeqNames(this, names)
    seqOrder <- typeOfSequenceOrdering(names, ...)
  }
  seqOrder
}, protected=TRUE)

setMethodS3("getSeqLengths", "SequenceContigsInterface", abstract=TRUE)

setMethodS3("hasSeqLengths", "SequenceContigsInterface", function(this, ...) {
  lens <- getSeqLengths(this, ...)
  !all(is.na(lens))
})

setMethodS3("getTotalSeqLength", "SequenceContigsInterface", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...)
  if (all(is.na(seqLengths))) return(NA_integer_)
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
  unames <- unique(names)
  map <- match(names, table=unames)
  cleanUnames <- sub(" .*", "", unames)
  names <- cleanUnames[map]
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

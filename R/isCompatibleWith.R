setMethodS3("isCompatibleWithBySeqNames", "default", function(this, other, mustWork=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':

  # Argument 'mustWork':
  mustWork <- Arguments$getLogical(mustWork)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  # Assert compatible sequence names
  verbose && enterf(verbose, "Checking whether %s is compatible with %s", class(this)[1L], class(other)[1L])

  verbose && print(verbose, this)
  verbose && print(verbose, other)

  if (inherits(this, "GenericDataFile")) {
    label <- getPathname(this);
  } else if (inherits(this, "GenericDataFileSet")) {
    label <- getPath(this);
  } else {
    label <- class(this)[1L]
  }

  if (inherits(other, "GenericDataFile")) {
    labelO <- getPathname(other);
  } else if (inherits(other, "GenericDataFileSet")) {
    labelO <- getPath(other);
  } else {
    labelO <- class(other)[1L]
  }

  verbose && enter(verbose, "Comparing sequence names")

  seqNames <- getSeqNames(this, unique=TRUE)
  verbose && printf(verbose, "%s: [%d] %s\n", class(this)[1L], length(seqNames), hpaste(sQuote(seqNames)))

  seqNamesO <- getSeqNames(other, unique=TRUE)
  verbose && printf(verbose, "%s: [%d] %s\n", class(other)[1L], length(seqNamesO), hpaste(sQuote(seqNamesO)))

  common <- intersect(seqNames, seqNamesO)
  verbose && printf(verbose, "In common: [%d] %s\n", length(common), hpaste(sQuote(common)))

  # Sanity check
  if (length(common) == 0L) {
    msg <- sprintf("The sequence names of the %s ('%s') and the %s set ('%s') are incompatible, because they have no names in common.", class(this)[1L], class(other)[1L], label, labelO);
    verbose && cat(verbose, msg)
    # Assertion?
    if (mustWork) throw(msg)
    verbose && exit(verbose)
    verbose && exit(verbose)
    return(FALSE)
  }

  verbose && exit(verbose)

  verbose && exit(verbose)

  TRUE
}, protected=TRUE) # isCompatibleWithBySeqNames()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SequenceContigsInterface
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("isCompatibleWithBySeqNames", "SequenceContigsInterface", function(this, other, unique=TRUE, mustWork=FALSE, ...) {
  res <- FALSE
  
  if (!inherits(other, "SequenceContigsInterface")) {
    msg <- sprintf("The 'other' object is not an SequenceContigsInterface class: %s", class(other)[1L])
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  names <- getSeqNames(this, unique=unique)
  namesO <- getSeqNames(other, unique=unique)

  ## Assert same names
  idxs <- match(names, namesO)
  nas <- is.na(idxs)
  if (any(nas)) {
    if (all(nas)) {
      msg <- "None of the sequence names matches."
    } else {
      unknown <- names[nas]
      msg <- sprintf("Some of the sequence names does not exist in target: [%d] %s", length(unknown), hpaste(sQuote(unknown)))
    }
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  ## Assert same order
  if (!all(diff(idxs) > 0, na.rm=TRUE)) {
    msg <- "The ordering of sequence names does not match."
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  TRUE
}, protected=TRUE)


setMethodS3("isCompatibleWithBySeqLengths", "SequenceContigsInterface", function(this, other, mustWork=FALSE, ...) {
  res <- FALSE
  
  if (!inherits(other, "SequenceContigsInterface")) {
    msg <- sprintf("The 'other' object is not an SequenceContigsInterface class: %s", class(other)[1L])
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  lens <- getSeqLengths(this)
  lensO <- getSeqLengths(other)

  ## Assert same lengths
  idxs <- match(lensO, lens)
  nas <- is.na(idxs)
  if (any(nas)) {
    if (all(nas)) {
      msg <- "None of the sequence lengths matches."
    } else {
      unknown <- lens[nas]
      msg <- sprintf("Some of the sequence lengths does not exist in target: [%d] %s", length(unknown), hpaste(sQuote(unknown)))
    }
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  ## Assert same order
  if (!all(diff(idxs) > 0, na.rm=TRUE)) {
    msg <- "The ordering of sequence lengths does not match."
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  TRUE
}, protected=TRUE)


setMethodS3("isCompatibleWithBySeqs", "SequenceContigsInterface", function(this, other, ...) {
  res <- isCompatibleWithBySeqNames(this, other, ...)
  if (!res) return(res)

  if (hasSeqLengths(this)) {
    res <- isCompatibleWithBySeqLengths(this, other, ...)
    if (!res) return(res)
  }

  TRUE
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File subclasses
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("isCompatibleWith", "FastaReferenceFile", function(this, other, ...) {
  isCompatibleWithBySeqNames(this, other, ...)
})

setMethodS3("isCompatibleWith", "AromaSeqDataFile", function(this, other, ...) {
  isCompatibleWithBySeqNames(this, other, ...)
})



setMethodS3("isCompatibleWith", "GcBaseFile", function(this, other, mustWork=FALSE, ...) {
  res <- NextMethod("isCompatibleWith")
  if (!isTRUE(res)) return(res)
  
  ## Assert same names
  names <- getSeqNames(other)
  idxs <- match(names, getSeqNames(this))
  nas <- is.na(idxs)
  if (any(nas)) {
    if (all(nas)) {
      msg <- "None of the sequence names matches."
    } else {
      unknown <- names[nas]
      msg <- sprintf("Some of the sequence names does not exist in target: [%d] %s", length(unknown), hpaste(sQuote(unknown)))
    }
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }

  ## Assert same order
  res <- all(diff(idxs) > 0, na.rm=TRUE)
  if (!res) {
    msg <- "The ordering of sequence names does not match."
    if (mustWork) throw(msg)
    attr(res, "reason") <- msg
    return(res)
  }
  
  TRUE
})

setMethodS3("isCompatibleWith", "Bowtie2IndexSet", function(this, other, mustWork=FALSE, ...) {
  if (isTopHat2IndexSet(this)) {
    if (inherits(other, "GtfDataFile")) {
      gtf <- other;
      fullname <- basename(getIndexPrefix(this));
      fullnameGTF <- getFullName(gtf);
      if (fullnameGTF != fullname) {
        msg <- sprintf("%s (%s) is incompatible with %s (%s) because their fullnames does not match: %s != %s", class(this)[1L], class(other)[1L], getPath(this), getPathname(gtf), fullname, fullnameGTF);
        if (mustWork) throw(msg);
        return(FALSE);
      }
    }
  } else {
    isCompatibleWithBySeqNames(this, other, mustWork=mustWork, ...)
  }
})

setMethodS3("isCompatibleWith", "BwaIndexSet", function(this, other, ...) {
  isCompatibleWithBySeqNames(this, other, ...)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# File set subclasses
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("isCompatibleWith", "AromaSeqDataFileSet", function(this, other, ...) {
  .stop_if_not(inherits(other, "SequenceContigsInterface"))

  for (ii in seq_along(this)) {
    df <- this[[ii]]
    res <- isCompatibleWith(df, other, ...)
    if (!res) return(res)
  }

  TRUE
})

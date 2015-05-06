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


setMethodS3("isCompatibleWith", "FastaReferenceFile", function(this, other, ...) {
  isCompatibleWithBySeqNames(this, other, ...)
})

setMethodS3("isCompatibleWith", "AromaSeqDataFile", function(this, other, ...) {
  isCompatibleWithBySeqNames(this, other, ...)
})

setMethodS3("isTopHat2IndexSet", "Bowtie2IndexSet", function(this, ...) {
  # AD HOC
  path <- getPath(this);
  grepl("tophat2", path, fixed=TRUE);
}, protected=TRUE)

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



############################################################################
# HISTORY:
# 2014-08-23
# o Added isCompatibleWith().
# o Added isCompatibleWithBySeqNames().
############################################################################


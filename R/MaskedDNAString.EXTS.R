setMethodS3("binTabulate", "MaskedDNAString", function(seq, bx, letters=c("A", "C", "G", "T"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'bx':
  if (!is.numeric(bx)) {
    stop("Argument 'bx' is not numeric: ", mode(bx));
  }
  if (any(is.infinite(bx))) {
    stop("Argument 'bx' must not contain Inf values.");
  }
  o <- order(bx);
  if (any(diff(o) != 1L)) {
    stop("Argument 'bx' is not ordered.");
  }

  # Argument 'letters':
  letters <- Arguments$getCharacters(letters);
  letters <- unlist(strsplit(letters, split=""), use.names=FALSE);


  # Pull out the actual sequence data
  # NB: Can this be done more efficient?
  # Alt 1: Not sure this is guaranteed to start with base #1. /HB 2012-10-17
  # x <- toString(seq);

  # Alt 2:
  seq <- unmasked(seq);
  # Get the offset of the sequences to be binned.  Should be zero
  # in most (all?) cases.  If non-zero, we'll adjust the binning
  # accordingly.
  offset <- seq@offset;
  .stop_if_not(offset >= 0L);

  # Adjust for sequence offset?
  if (offset != 0L) {
    # Adjusting the bins, rather than the sequence offset,
    # should be more efficient since they are fewer and it
    # has only to be once.
    bx <- bx - offset;
  }

  # Coerce to a character string.
  x <- toString(seq);

  # Not needed anymore
  rm(seq);

  # Coerce into raw
  x <- charToRaw(x);

  # Tabulate
  counts <- matrix(0L, nrow=length(bx)-1L, ncol=length(letters));
  colnames(counts) <- letters;
  for (letter in letters) {
    z <- which(x == charToRaw(letter));
    counts[,letter] <- binCounts(z, bx=bx, ...);
    rm(z);
  }

  # Not needed anymore
  rm(x);

  counts;
}) # binTabulate()

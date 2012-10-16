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
  x <- toString(seq);
  x <- charToRaw(x);

  # Tabulate
  counts <- matrix(0L, nrow=length(bx)-1L, ncol=length(letters));
  colnames(counts) <- letters;
  for (letter in letters) {
    z <- which(x == charToRaw(letter));
    counts[,letter] <- binCounts(z, bx=bx, ...);
    rm(z);
  }
  rm(x);

  counts;
}) # binTabulate()


############################################################################
# HISTORY:
# 2012-10-16
# o Added binTabulate() for MaskedDNAString.
# o Created.
############################################################################

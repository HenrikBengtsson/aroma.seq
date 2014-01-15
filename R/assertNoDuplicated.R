assertNoDuplicated <- function(x, .name=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument '.name':
  if (is.null(.name)) {
    .name <- as.character(deparse(substitute(x)));
  }

  # Argument 'x':
  x <- Arguments$getVector(x, .name=.name);

  # Nothing todo?
  if (length(x) == 1L) return(x);


  # No duplicates?
  if (!anyDuplicated(x)) return(x);

  # Has duplicates!
  dups <- x[duplicated(x)];
  throw(sprintf("Detected %d duplicated entries in argument '%s': %s", length(dups), .name, hpaste(sQuote(dups))));
} # assertNoDuplicated()


############################################################################
# HISTORY:
# 2014-01-14
# o Created.
############################################################################

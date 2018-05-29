hasCommas <- function(s, ...) {
  (regexpr(",", s, fixed=TRUE) != -1L)
}

hasSpaces <- function(s, ...) {
  (regexpr(" ", s, fixed=TRUE) != -1L)
}

assertNoCommas <- function(s, ...) {
  .stop_if_not(!any(hasCommas(s)))
}

assertNoSpaces <- function(s, ...) {
  .stop_if_not(!any(hasSpaces(s)))
}

assertNoDuplicated <- function(x, .name=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument '.name':
  if (is.null(.name)) {
    .name <- as.character(deparse(substitute(x)))
  }

  # Argument 'x':
  x <- Arguments$getVector(x, .name=.name)

  # Nothing todo?
  if (length(x) == 1L) return(x)


  # No duplicates?
  if (!anyDuplicated(x)) return(x)

  # Has duplicates!
  dups <- x[duplicated(x)]
  throw(sprintf("Detected %d duplicated entries in argument '%s': %s", length(dups), .name, hpaste(sQuote(dups))))
} # assertNoDuplicated()


setMethodS3("getTuxedoOption", "Arguments", function(static, value, ..., .name=NULL) {
  .name <- as.character(deparse(substitute(value)))
  nok <- hasCommas(value)
  if (any(nok)) {
    throw(sprintf("Tuxedo (Bowtie, TopHat, ...) option %s must not contain commas: %s", sQuote(.name), hpaste(sQuote(value[nok]), collapse=", ")))
  }
  nok <- hasSpaces(value)
  if (any(nok)) {
    throw(sprintf("Tuxedo (Bowtie, TopHat, ...) option %s must not contain spaces: %s", sQuote(.name), hpaste(sQuote(value[nok]), collapse=", ")))
  }
  value
}, static=TRUE)

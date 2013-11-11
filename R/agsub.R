agsub <- function(pattern, replacement, x, ..., default=NA_character_, as=c("matrix", "list", "data.frame")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pattern'
  if (length(pattern) != 1L) {
    throw("Argument 'pattern' must be a single character string.");
  }

  # Argument 'replacement'
  if (!is.vector(replacement)) {
    throw("Argument 'replacement' must be a vector.");
  }

  # Argument 'default'
  if (length(default) != 1L) {
    throw("Argument 'default' must be a single character string.");
  }

  # Argument 'as'
  as <- match.arg(as);


  # Coerce vector to list
  if (!is.list(replacement)) {
    replacement <- as.list(replacement);
  }

  # Expand default
  default <- rep(default, length=length(x));

  # Allocate result
  res <- vector("list", length=length(replacement));
  names(res) <- names(replacement);

  # Identifying matching strings, iff any
  idxs <- which(regexpr(pattern, x, ...) != -1L);
  match <- (length(idxs) > 0L);

  # Assign according to replacement vector
  for (kk in seq_along(replacement)) {
    resKK <- default;
    if (match) {
      resKK[idxs] <- gsub(pattern, replacement[[kk]], x[idxs], ...);
    }
    res[[kk]] <- resKK;
  } # for (kk ...)

  # Coerce results?
  if (as == "matrix") {
    res <- unlist(res, use.names=FALSE);
    dim(res) <- c(length(x), length(replacement));
    colnames(res) <- names(replacement);
    rownames(res) <- names(x);
  } else if (as == "data.frame") {
    res <- as.data.frame(res, check.names=FALSE, stringsAsFactors=FALSE);
    colnames(res) <- names(replacement);
#    rownames(res) <- names(x);
  }

  res;
} # agsub()


############################################################################
# HISTORY:
# 2013-11-10
# o Added agsub(). Will probably end up in 'R.utils' at some point.
# o Created.
############################################################################

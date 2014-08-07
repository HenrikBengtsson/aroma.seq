findHTSeq <- function(..., command=c("htseq-count", "htseq-qa")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command);

  # Aroma-specific variable
  path <- getExternalHome("HTSEQ_HOME");

  versionPattern <- c(".*version ([0-9.]+(|p[0-9]+)).*");
  res <- findExternal(command=command, path=path, versionPattern=versionPattern, ...);

  # Update version format '0.5.4p3' to '0.5.4-3'
  if (!is.null(res)) {
    ver <- attr(res, "version");
    ver <- gsub("p", "-", ver, fixed=TRUE);
    attr(res, "version") <- ver;
  }

  res;
} # findHTSeq()

############################################################################
# HISTORY:
# 2014-03-09 [HB]
# o Now findHTSeq() returns versions in format '0.5.4-3' not '0.5.4p3',
# 2014-01-24 [HB]
# o BUG FIX: findHTSeq() failed to identify the version.
# 2013-06-20 [HB]
# o Renamed from findHtseq() to findHTSeq().
# 2013-05-31 [TT]
# o Created.
############################################################################
